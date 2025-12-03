import os
import re
import sys
import json
import math
import shutil
import argparse
import warnings
import numpy as np
from pyatb.easy_use.stru_analyzer import read_abacus_stru
from ase import Atoms
from typing import List, Tuple

def _parse_upf_valence(upf_path):
    """从 UPF 文件解析 z_valence（优先 v2 的 PP_HEADER，其次兼容部分 v1 写法）。"""
    with open(upf_path, "r", encoding="utf-8", errors="ignore") as f:
        head = f.read(200000)

    m = re.search(r'z[_\s]?valence\s*=\s*"?\s*([0-9]+(?:\.[0-9]*)?)', head, re.IGNORECASE)  # v2
    if m: return float(m.group(1))

    m = re.search(r'<\s*z[_\s]?valence\s*>\s*([0-9]+(?:\.[0-9]*)?)\s*<\s*/\s*z[_\s]?valence\s*>',
                  head, re.IGNORECASE)  # v1
    if m: return float(m.group(1))

    m = re.search(r'^\s*z[_\s]?valence\s+([0-9]+(?:\.[0-9]*)?)\s*$', head,
                  re.IGNORECASE | re.MULTILINE)  # 行式
    if m: return float(m.group(1))

    m = re.search(r'valence\s+charge[^0-9]*([0-9]+(?:\.[0-9]*)?)', head, re.IGNORECASE)  # 兜底
    if m: return float(m.group(1))

    raise ValueError(f"无法从 UPF 中解析 z_valence: {upf_path}")


def _first_occurrence_species(ase_stru):
    """按结构中首次出现顺序生成元素列表。"""
    order, seen = [], set()
    for s in ase_stru.get_chemical_symbols():
        if s not in seen:
            order.append(s)
            seen.add(s)
    return order

def parse_input(input_text):
    variables = {}
    for line in input_text.split("\n"):
        if line.strip() == "" or line.strip().startswith("#"):
            continue
        line_content = line.split("#")[0].strip()
        if line_content and " " in line_content:
            parts = line_content.split(maxsplit=1)
            if len(parts) == 2:
                key, value = parts
                variables[key] = value
    return variables

def parse_input_file(file_path):
    try:
        with open(file_path, 'r') as file:
            input_text = file.read()
            return parse_input(input_text)
    except FileNotFoundError:
        print(f"文件 {file_path} 未找到。")
        return {}

def extract_properties_from_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    max_val = data.get("max_val", None)
    e_fermi = data.get("fermi_energy", None)
    band_gap = data.get("band_gap", None)
    n_elec = int(data.get("num_electrons", 0))
    
    # n_occu = 一半取整
    n_occu = math.floor(n_elec / 2)

    return max_val, e_fermi, band_gap, n_elec, n_occu



def kpath_generator(ase_stru: Atoms, kline_density = 0.01, dim = '3', tolerance=5e-4, knum=0, kpath = None):
    """
    Generate k-points along a high-symmetry path in reciprocal space.

    Args:
        ase_stru (Atoms): The atomic structure.
        kline_density (float, optional): The density of k-points along the path. Defaults to 0.01.
        dim (str, optional): The dimension of the system. Can be '2', '3', or a string of three integers separated by spaces. Defaults to '3'.
        kpath (list, optional): The path in reciprocal space. Defaults to None.
        tolerance (float, optional): The tolerance for determining the path. Defaults to None.

    Returns:
        tuple: A tuple containing three lists: kpt_output, kpoint_label, and kpoint_num_in_line.
            - kpt_output (list): The generated k-points with their corresponding densities and labels.
            - kpoint_label (list): The labels of the k-points.
            - kpoint_num_in_line (list): The number of k-points in each line.

    """
    # 根据维度 dim 设置路径参数 pbc
    if dim == '2':
        bandpath = ase_stru.cell.bandpath(path = kpath, density=kline_density , eps = tolerance, pbc=[1, 1, 0])
    elif dim == '3':
        bandpath = ase_stru.cell.bandpath(path = kpath, density=kline_density , eps = tolerance)
    elif dim == 'x':
        bandpath = ase_stru.cell.bandpath(path = kpath, density=kline_density , eps = tolerance, pbc=[1, 0, 0])
    elif dim == 'y':
        bandpath = ase_stru.cell.bandpath(path = kpath, density=kline_density , eps = tolerance, pbc=[0, 1, 0])
    else :
        if ' ' in dim:
            dim = [int(dimx) for dimx in dim.split()]
            bandpath = ase_stru.cell.bandpath(path = kpath, density=kline_density , eps = tolerance, pbc=dim)
        else:
            return

    path_label = bandpath.path
    # 使用正则表达式匹配字母和数字组合，例如 'X1'
    pattern = re.findall(r'[A-Z][0-9]*', path_label)
    path_label_array = list(pattern)
    cleaned_labels = [label for label in pattern if label != ',']

    special_points = bandpath.special_points
    # 对每个坐标进行检查和变换 只对不在 [-1, 1] 范围内的坐标应用 mod 1 操作
    shifted_points = {
        key: value if np.all((value >= -1) & (value <= 1)) else np.mod(value, 1)
        for key, value in special_points.items()
    }
    special_points = shifted_points

    kpt_output = []
    kpoint_label = []
    kpoint_num_in_line = []

    rec_lat_cell = ase_stru.cell.reciprocal()
    rec_lat_matrix  = rec_lat_cell[:]
    # print(rec_lat_matrix)
    for i in range(len(path_label_array)):
        # 在这里使用 i 和 path_label_array[i] 进行操作
        label = path_label_array[i]
        label_next = path_label_array[i+1] if i+1 < len(path_label_array) else None
        coordinates = special_points.get(label)
        coordinates_next = special_points.get(label_next) if label_next else None  # 如果获取不到值也置为None
        if coordinates is not None and coordinates_next is not None:
            if knum == 0: # 如果knum为0，则计算两点之间的k点数
                # 计算两点之间的距离
                k_real = coordinates @ rec_lat_matrix
                k_real_next = coordinates_next @ rec_lat_matrix
                distance_in_reciprocal = np.linalg.norm(k_real - k_real_next)
                distance = distance_in_reciprocal
                # print(f"Distance in reciprocal space is {distance_in_reciprocal}")
                # 计算两点之间的k点数，如果小于3则置为3
                density = max(int(distance * (2* np.pi) / kline_density ), 3)
                kpt_output.append(f"{'  '.join([f'{coord: .10f}' for coord in coordinates])}  {format(density, '<4')}   # {label}")
                kpoint_label.append(f"{label}   ")
                kpoint_num_in_line.append(f"{density}  ")
            else: # 如果knum不为0，则直接使用knum设置k点数
                kpt_output.append(f"{'  '.join([f'{coord: .10f}' for coord in coordinates])}  {format(knum, '<4')}   # {label}")
                kpoint_label.append(f"{label}   ")
                kpoint_num_in_line.append(f"{knum}  ")
        elif coordinates is None:
            # print("no coord this line")
            pass
        else:
            kpt_output.append(f"{'  '.join([f'{coord: .10f}' for coord in coordinates])}  {format('1', '<4')}   # {label}")
            kpoint_label.append(f"{label}   ")
            kpoint_num_in_line.append(f"{ format('1', '<4') }")
    
    kpt_output.append(f' kpoint_label{" " * 20}{",".join(cleaned_labels)}')
    # kpt_output.append(f" # Auto-generated K-path: [{path_label}]")

    return kpt_output, kpoint_label, kpoint_num_in_line


def generate_input_init(nspin, e_fermi, out_suffix_path, lattice_constant, lattice_vectors, max_kpoint_num):
    """
    Generate the input text for initializing a simulation.

    Args:
        nspin (int): Number of spins.
        e_fermi (float): Fermi energy.
        out_suffix_path (str): Output suffix path.
        lattice_constant (float): Lattice constant.
        lattice_vectors (list): List of lattice vectors.
        max_kpoint_num (int): Maximum number of k-points.

    Returns:
        str: The generated input text.
    """
    if nspin == 2:
        HR_file = str(os.path.join(out_suffix_path, "data-HR-sparse_SPIN0.csr")) + ", " + str(os.path.join(out_suffix_path, "data-HR-sparse_SPIN1.csr"))
    else:
        HR_file = str(os.path.join(out_suffix_path, "data-HR-sparse_SPIN0.csr"))
    
    input_text = f"INPUT_PARAMETERS\n{{\n"
    input_text += f"    nspin                           {nspin}\n"
    input_text += f"    package                         ABACUS\n"
    input_text += f"    fermi_energy                    {e_fermi}\n"
    input_text += f"    fermi_energy_unit               eV\n"
    input_text += f"    HR_route                        {HR_file}\n"
    input_text += f"    SR_route                        {out_suffix_path}/data-SR-sparse_SPIN0.csr\n"
    input_text += f"    rR_route                        {out_suffix_path}/data-rR-sparse.csr\n"
    input_text += f"    HR_unit                         Ry\n"
    input_text += f"    rR_unit                         Bohr\n"
    input_text += f"    max_kpoint_num                  {max_kpoint_num}\n"
    input_text += f"}}\n\n"
    input_text += f"LATTICE\n{{\n"
    input_text += f"    lattice_constant                {lattice_constant}\n"
    input_text += f"    lattice_constant_unit           Angstrom\n"
    input_text += f"    lattice_vector\n"
    input_text += f"    {lattice_vectors[0][0]}  {lattice_vectors[0][1]}  {lattice_vectors[0][2]}\n"
    input_text += f"    {lattice_vectors[1][0]}  {lattice_vectors[1][1]}  {lattice_vectors[1][2]}\n"
    input_text += f"    {lattice_vectors[2][0]}  {lattice_vectors[2][1]}  {lattice_vectors[2][2]}\n"
    input_text += f"}}\n"
    return input_text


def generate_hamgnn_input_init(nspin, e_fermi, out_suffix_path, lattice_constant, lattice_vectors, max_kpoint_num):
    """
    Generate the input text for initializing a simulation.

    Args:
        nspin (int): Number of spins.
        e_fermi (float): Fermi energy.
        out_suffix_path (str): Output suffix path.
        lattice_constant (float): Lattice constant.
        lattice_vectors (list): List of lattice vectors.
        max_kpoint_num (int): Maximum number of k-points.

    Returns:
        str: The generated input text.
    """
    if nspin == 2:
        HR_file = str(os.path.join(out_suffix_path, "H0.csr")) + ", " + str(os.path.join(out_suffix_path, "H1.csr"))
    else:
        HR_file = str(os.path.join(out_suffix_path, "H.csr"))
    
    input_text = f"INPUT_PARAMETERS\n{{\n"
    input_text += f"    nspin                           {nspin}\n"
    input_text += f"    package                         ABACUS\n"
    input_text += f"    fermi_energy                    {e_fermi}\n"
    input_text += f"    fermi_energy_unit               eV\n"
    input_text += f"    HR_route                        {HR_file}\n"
    input_text += f"    SR_route                        {out_suffix_path}/S.csr\n"
    input_text += f"    rR_route                        {out_suffix_path}/R.csr\n"
    input_text += f"    HR_unit                         eV\n"
    input_text += f"    rR_unit                         Bohr\n"
    input_text += f"    max_kpoint_num                  {max_kpoint_num}\n"
    input_text += f"}}\n\n"
    input_text += f"LATTICE\n{{\n"
    input_text += f"    lattice_constant                {lattice_constant}\n"
    input_text += f"    lattice_constant_unit           Angstrom\n"
    input_text += f"    lattice_vector\n"
    input_text += f"    {lattice_vectors[0][0]}  {lattice_vectors[0][1]}  {lattice_vectors[0][2]}\n"
    input_text += f"    {lattice_vectors[1][0]}  {lattice_vectors[1][1]}  {lattice_vectors[1][2]}\n"
    input_text += f"    {lattice_vectors[2][0]}  {lattice_vectors[2][1]}  {lattice_vectors[2][2]}\n"
    input_text += f"}}\n"
    return input_text


def generate_input_band(input_text, ase_stru: Atoms, kline_density=100, dim=3, kmode='line', tolerance=5e-4, knum=0, kpath = None):
    """
    Generates the input for band structure calculation.

    Args:
        input_text (str): The existing input text.
        ase_stru (Atoms): The ASE structure object.
        kline_density (int, optional): The density of k-points along the band line. Defaults to 100.
        dim (int, optional): The dimension of the k-point grid. Defaults to 3.
        kmode (str, optional): The k-point mode. Defaults to 'line'.
        kpath (list, optional): The list of high symmetry k-points. Defaults to None.
        tolerance (float, optional): The tolerance for determining high symmetry k-points. Defaults to 5e-4.

    Returns:
        str: The updated input text.
    """
    lattice_vectors = ase_stru.get_cell()
    # 获取band_line的行数
    input_text += f"\nBAND_STRUCTURE\n"
    input_text += f"{{\n"
    input_text += f"    wf_collect                      0\n"
    if kmode in ["mp", "mesh"]:
        integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.03)  # band mesh 默认密度
        input_text += f"    kpoint_mode               mp\n"
        input_text += f"    mp_grid                   {' '.join(map(str, integrate_grid))}\n"
    else:
        band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
        num_lines = len(band_line) - 1 
        input_text += f"    kpoint_mode               {kmode}\n"
        input_text += f"    kpoint_num                {num_lines}\n"
        input_text += f"    high_symmetry_kpoint\n"
        for line in band_line:
            input_text += f"   {line}\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text

def generate_input_fatband(input_text, ase_stru: Atoms, kline_density=100, dim=3, kmode='line',  tolerance=5e-4,  knum=0, kpath = None,  band_range=None, allbands=None):
    """
    Generates input for the FAT_BAND calculation in pyatb.

    Args:
        input_text (str): The existing input text.
        ase_stru (Atoms): The ASE structure object.
        kline_density (int, optional): The density of k-points along the k-path. Defaults to 100.
        dim (int, optional): The dimension of the system. Defaults to 3.
        kpath (list, optional): The k-path. Defaults to None.
        tolerance (float, optional): The tolerance for symmetry analysis. Defaults to None.
        band_range (list, optional): The range of bands to calculate. Defaults to None.
        allbands (bool, optional): Whether to calculate all bands. Defaults to None.

    Returns:
        str: The updated input text.
    """
    band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)

    lattice_vectors = ase_stru.get_cell()

    # 获取band_line的行数
    num_lines = len(band_line) - 1 
    input_text += f"\nFAT_BAND\n"
    input_text += f"{{\n"
    input_text += f"    band_range                      {band_range[0]}  {band_range[1]}\n"
    input_text += f"    stru_file                       STRU\n"

    if kmode in ["mp", "mesh"]:
        integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.015)  # band mesh 默认密度
        input_text += f"    kpoint_mode               mp\n"
        input_text += f"    mp_grid                   {' '.join(map(str, integrate_grid))}\n"
    else:
        band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
        num_lines = len(band_line) - 1 
        input_text += f"    kpoint_mode               {kmode}\n"
        input_text += f"    kpoint_num                {num_lines}\n"
        input_text += f"    high_symmetry_kpoint\n"
        for line in band_line:
            input_text += f"   {line}\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text


def generate_input_spintexture(input_text, ase_stru: Atoms, kline_density=100, dim=3, kmode='line', tolerance=5e-4,  knum=0, kpath = None,  band_range=None, allbands=None):
    band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
    lattice_vectors = ase_stru.get_cell()

    # 获取band_line的行数
    num_lines = len(band_line) - 1 
    input_text += f"\nSPIN_TEXTURE\n"
    input_text += f"{{\n"
    input_text += f"    band_range                      {band_range[0]}  {band_range[1]}\n"

    if kmode in ["mp", "mesh"]:
        integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.012)  # band mesh 默认密度
        input_text += f"    kpoint_mode               mp\n"
        input_text += f"    mp_grid                   {' '.join(map(str, integrate_grid))}\n"
    else:
        band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
        num_lines = len(band_line) - 1
        input_text += f"    kpoint_mode               {kmode}\n"
        input_text += f"    kpoint_num                {num_lines}\n"
        input_text += f"    high_symmetry_kpoint\n"
        for line in band_line:
            input_text += f"   {line}\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text

def generate_input_bandunfold(input_text, ase_stru: Atoms, kline_density=100, dim=3, kmode='line', tolerance=5e-4,  knum=0, m_matrix='1 0 0   0 1 0   0 0 1', kpath = None,  band_range=None, allbands=None):
    """
    生成用于能带反折叠计算的输入。

    参数:
        input_text (str): 现有的输入文本。
        ase_stru (Atoms): ASE 结构对象。
        kline_density (int, optional): 沿 k 路径的 k 点密度。默认为 100。
        dim (int, optional): 系统的维度。默认为 3。
        kpath (list, optional): 沿 k 路径的 k 点列表。默认为 None。
        tolerance (float, optional): 对称性分析的容差。默认为 None。
        band_range (str, optional): 要展开的能带的范围。默认为 None。
        allbands (int, optional): 总带数。默认为 None。

    返回:
        str: 更新后的输入文本。
    """

    lattice_vectors = ase_stru.get_cell()
    band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
    # 获取band_line的行数
    num_lines = len(band_line) - 1 
    input_text += f"\nBANDUNFOLDING\n"
    input_text += f"{{\n"
    input_text += f"    band_range                      {band_range[0]}  {band_range[1]}\n"
    input_text += f"    ecut                            10\n"
    input_text += f"    stru_file                       STRU\n"
    input_text += f"    m_matrix                        {m_matrix}\n"
    if kmode in ["mp", "mesh"]:
        integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.012)  # band mesh 默认密度
        input_text += f"    kpoint_mode               mp\n"
        input_text += f"    mp_grid                   {' '.join(map(str, integrate_grid))}\n"
    else:
        band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
        num_lines = len(band_line) - 1
        input_text += f"    kpoint_mode               {kmode}\n"
        input_text += f"    kpoint_num                {num_lines}\n"
        input_text += f"    high_symmetry_kpoint\n"
        for line in band_line:
            input_text += f"   {line}\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text

def generate_input_pdos(input_text, dim, lattice_vectors, e_fermi, energy_range):
    """
    根据提供的参数生成输入文本来计算PDOS。

    :param input_text: 基础的输入文本模板。
    :param e_fermi: 费米能级。
    :param energy_range: 能量范围参数，可以是(e_low, e_high)的形式，也可以是单一值代表与费米能级的距离。
    :return: 修改后包含能量范围的输入文本。
    """
    e_low, e_high = e_fermi + energy_range[0], e_fermi + energy_range[1]

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.10)  # pdos 默认密度

    input_text += f"\nPDOS\n"
    input_text += f"{{\n"
    input_text += f"    stru_file      STRU\n"
    input_text += f"    e_range        {e_low} {e_high}\n"
    input_text += f"    de             0.01\n"
    input_text += f"    sigma          0.10\n"
    input_text += f"    kpoint_mode    mp\n"
    input_text += f"    mp_grid        {' '.join(map(str, integrate_grid))}\n"
    input_text += f"}}\n"

    return input_text

def generate_input_findnodes(input_text, dim, lattice_vectors, e_fermi, energy_range):

    e_low, e_high = e_fermi + energy_range[0], e_fermi + energy_range[1]

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.06)  # FIND NODES 默认密度

    input_text += f"\nFIND_NODES\n"
    input_text += f"{{\n"
    input_text += f"    energy_range    {e_low} {e_high}\n"
    input_text += f"    k_start  0.0  0.0 -0.2\n"
    input_text += f"    k_vect1  0.0  0.0  0.0\n"
    input_text += f"    k_vect2  0.0  0.0  0.0\n"
    input_text += f"    k_vect2  0.0  0.0  0.4\n"
    input_text += f"    initial_grid              {' '.join(map(str, integrate_grid))}\n"
    input_text += f"    initial_threshold         0.01\n"
    input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
    input_text += f"    adaptive_threshold        0.001\n"
    input_text += f"}}\n"

    return input_text


def generate_input_optical(input_text, dim, lattice_vectors, n_occu, omega_range, mp_density):
    """
    生成光学输入的函数。

    Args:
        input_text (str): 输入文本。
        dim (int): 维度。
        lattice_vectors (list): 晶格矢量。
        n_occu (int): 占据态数目。
        *omega_range (float): omega范围参数。

    Returns:
        str: 生成的光学输入文本。
    """
    if not omega_range:
        omega_low = 0
        omega_high = 10
    elif len(omega_range) == 1:
        omega_low = 0
        omega_high = omega_range[0]
    elif len(omega_range) == 2:
        omega_low, omega_high = omega_range
    else:
        raise ValueError("omega范围参数应该是单一值或一对值。")

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.05)  # optical 默认密度

    input_text += f"\nOPTICAL_CONDUCTIVITY\n"
    input_text += f"{{\n"
    input_text += f"    occ_band       {n_occu}\n"
    input_text += f"    omega          {omega_low} {omega_high}\n"
    input_text += f"    domega         0.01\n"
    input_text += f"    eta            0.10\n"
    input_text += f"    grid           {' '.join(map(str, integrate_grid))}\n"
    input_text += f"}}\n"

    return input_text

def generate_input_jdos(input_text, dim, lattice_vectors, n_occu, omega_range, mp_density):
    if not omega_range:
        omega_low = 0
        omega_high = 5
    elif len(omega_range) == 1:
        omega_low = 0
        omega_high = omega_range[0]
    elif len(omega_range) == 2:
        omega_low, omega_high = omega_range
    else:
        raise ValueError("omega范围参数应该是单一值或一对值。")

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, mp_density)  # JDOS 默认密度

    input_text += f"\nJDOS\n"
    input_text += f"{{\n"
    input_text += f"    occ_band       {n_occu}\n"
    input_text += f"    omega          {omega_low} {omega_high}\n"
    input_text += f"    domega         0.01\n"
    input_text += f"    eta            0.10\n"
    input_text += f"    grid           {' '.join(map(str, integrate_grid))}\n"
    input_text += f"}}\n"

    return input_text

def generate_input_shift(input_text, dim, lattice_vectors, n_occu, omega_range, mp_density):
    if not omega_range:
        omega_low = 0
        omega_high = 5
    elif len(omega_range) == 1:
        omega_low = 0
        omega_high = omega_range[0]
    elif len(omega_range) == 2:
        omega_low, omega_high = omega_range
    else:
        raise ValueError("omega范围参数应该是单一值或一对值。")

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, mp_density)  # SHIFT CURRENT 默认密度

    input_text += f"\nSHIFT_CURRENT\n"
    input_text += f"{{\n"
    input_text += f"    occ_band       {n_occu}\n"
    input_text += f"    omega          {omega_low} {omega_high}\n"
    input_text += f"    domega         0.01\n"
    input_text += f"    eta            0.10\n"
    input_text += f"    grid           {' '.join(map(str, integrate_grid))}\n"
    input_text += f"}}\n"

    return input_text


def get_k_mesh_from_dim(lattice_vectors, dim, mesh_density):
    """
    Calculate the k-mesh parameters based on the lattice vectors, dimension, and mesh density.

    Parameters:
    lattice_vectors (array-like): The lattice vectors of the crystal.
    dim (int): The dimension of the crystal (2 or 3).
    mesh_density (float): The desired density of the k-mesh.

    Returns:
    tuple: A tuple containing the integrate grid, adaptive grid, and periodic boundary conditions (pbc).
    integrate_grid (list): The number of grid points in each direction for integration.
    adaptive_grid (list): The number of grid points in each direction for adaptive sampling.
    pbc (list): The periodic boundary conditions for each direction.

    Raises:
    ValueError: If the input values are invalid.
    """
    adaptive_density_ratio = 10
    # 解析 dim 并设置 pbc
    if dim == '3':
        pbc = [1, 1, 1]
    elif dim == '2':
        pbc = [1, 1, 0]
    elif dim == 'x':
        pbc = [1, 0, 0]
    elif dim == 'y':
        pbc = [0, 1, 0]
    elif len(dim) == 3 and all(x in ['0', '1'] for x in dim):
        pbc = [int(x) for x in dim]
    else:
        raise ValueError("Invalid value for dim. Please provide 3/2/x/y/z or a 3 string to define the pbc, e.g. 101 to define the periodic condition along x and z.")

    # 检查 pbc 的长度是否为 3
    if len(pbc) != 3:
        raise ValueError("Invalid length for pbc, must be 3")

    # 计算每个方向上的格点数
    reciprocal_lattice_vectors = 2 * np.pi * np.linalg.inv(lattice_vectors).T
    lengths = np.linalg.norm(reciprocal_lattice_vectors, axis=1)
    
    k_mesh = [max(1, int(length / mesh_density)) if p == 1 else 1 for length, p in zip(lengths, pbc)]

    # 根据 pbc 中 1 的数量设置网格参数
    if sum(pbc) == 3:
        integrate_grid = k_mesh
        adaptive_grid = [max(1, int(k / adaptive_density_ratio)) for k in k_mesh]
    elif sum(pbc) == 2:
        integrate_grid = [k if p == 1 else 1 for k, p in zip(k_mesh, pbc)]
        adaptive_grid = [max(1, int(k / adaptive_density_ratio)) if p == 1 else 1 for k, p in zip(k_mesh, pbc)]
    elif sum(pbc) == 1:
        integrate_grid = [k if p == 1 else 1 for k, p in zip(k_mesh, pbc)]
        adaptive_grid = [max(1, int(k / adaptive_density_ratio)) if p == 1 else 1 for k, p in zip(k_mesh, pbc)]
    else:
        raise ValueError("Invalid pbc configuration, at least one dimension should be periodic")

    # 调整 adaptive_grid，确保最小值不小于 5，最大值不超过 20
    adaptive_grid = [min(max(k, 5), 20) if p == 1 else 1 for k, p in zip(adaptive_grid, pbc)]
    
    return integrate_grid, adaptive_grid, pbc



def generate_input_bcd(input_text, dim, lattice_vectors, e_fermi, energy_range):

    e_low, e_high = e_fermi + energy_range[0], e_fermi + energy_range[1]

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.012)  # bcd 密度
    input_text += f"\nBERRY_CURVATURE_DIPOLE\n"
    input_text += f"{{\n"
    input_text += f"    omega              {e_low} {e_high}\n"
    input_text += f"    domega             0.001\n"
    input_text += f"    grid               {' '.join(map(str, integrate_grid))}\n"
    # input_text += f"    integrate_mode            Grid\n"
    # input_text += f"    integrate_grid            {' '.join(map(str, integrate_grid))}\n"
    # input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
    # input_text += f"    adaptive_grid_threshold   2000\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text

def generate_input_cpge(input_text, dim, lattice_vectors, e_fermi, energy_range):
    # 检查energy_range的长度，来决定是直接使用还是基于费米能级计算

    e_low, e_high = e_fermi + energy_range[0], e_fermi + energy_range[1]

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.01)  # cpge 密度

    input_text += f"\nCPGE\n"
    input_text += f"{{\n"
    input_text += f"    omega          {e_low} {e_high}\n"
    input_text += f"    domega         0.001\n"
    input_text += f"    integrate_mode            Grid\n"
    input_text += f"    integrate_grid            {' '.join(map(str, integrate_grid))}\n"
    input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
    input_text += f"    adaptive_grid_threshold   2000\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text


def generate_input_ahc(input_text, dim, lattice_vectors, mp_density=0.008):

    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, mp_density)  # AHC 默认密度 0.008


    input_text += f"\nAHC\n"
    input_text += f"{{\n"
    input_text += f"    integrate_mode            Grid\n"
    input_text += f"    integrate_grid            {' '.join(map(str, integrate_grid))}\n"
    input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
    input_text += f"    adaptive_grid_threshold   1000\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text

def generate_input_polarization(input_text, dim, lattice_vectors, n_occu,
                                ase_stru, stru_path, valence_str, mp_density,
                                stru_file_field="STRU"):
    """
    生成极化输入（POLARIZATION）——仿写 generate_input_optical 的风格与拼接方式。

    Args:
        input_text (str): 输入文本（将把 POLARIZATION 区块追加其后）。
        dim (int or str): 维度/周期性定义，透传给 get_k_mesh_from_dim（如 '3'/'2'/'101' 等）。
        lattice_vectors (array-like): 晶格矢量。
        n_occu (int): 占据态数目（occ_band）。
        ase_stru: ASE Atoms 对象（需包含 info['pp'] 映射）。
        stru_path (str): 对应 STRU 文件路径（用于解析相对 PP 路径基准）。
        mp_density (float): k 网格密度。
        stru_file_field (str): POLARIZATION 中的 stru_file 字段（默认 "STRU"）。
        valence_str (str): "auto"（默认自动解析）或形如 "13 3 2" 的手动字符串。
    生成极化输入（POLARIZATION）——仿写 generate_input_optical 的风格与拼接方式。
    - 当 valence_str=="auto" 且无法从 ase_stru.info['pp'] 解析 valence 时，给出警告并继续运行，
      此时 valence_line 设为占位字符串 "set the valence for the element in your STRU"。

    Returns:
        str: 追加了 POLARIZATION 区块后的文本。

    """
    # 计算 k 网格（与你的 optical 风格保持一致）
    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(np.array(lattice_vectors), dim, mp_density)
    nk1, nk2, nk3 = map(int, integrate_grid)

    # 生成 valence_e 与 atom_type
    use_placeholder = False
    placeholder_valence_line = "  #  number of valence electrons for the elements; it should match the Pseudopotentials, for example: 12 6"

    if isinstance(valence_str, str) and valence_str.strip().lower() != "auto":
        # 手动模式：按空格切分并转为整数
        tokens = valence_str.strip().split()
        if not tokens:
            warnings.warn("valence_str 提供了空字符串，将使用占位提示。")
            use_placeholder = True
            # atom_type 仍可依据结构推断
            species_order = _first_occurrence_species(ase_stru)
            atom_type = len(species_order)
            valence_line = placeholder_valence_line
        else:
            try:
                valences = [int(round(float(t))) for t in tokens]
                atom_type = len(valences)
                valence_line = ' '.join(map(str, valences))
            except ValueError:
                warnings.warn(f"valence_str 解析失败（应为空格分隔的数字）：{valence_str}；将使用占位提示。")
                use_placeholder = True
                species_order = _first_occurrence_species(ase_stru)
                atom_type = len(species_order)
                valence_line = placeholder_valence_line
    else:
        # 自动模式：尝试从 ase_stru.info['pp'] 解析
        species_order = _first_occurrence_species(ase_stru)
        atom_type = len(species_order)

        pp_map = getattr(ase_stru, "info", {}).get("pp", None)
        if not isinstance(pp_map, dict) or not pp_map:
            warnings.warn("未找到有效的 ase_stru.info['pp']（缺失/不是字典/为空）。将使用占位提示继续运行。")
            use_placeholder = True
            valence_line = placeholder_valence_line
        else:
            base_dir = os.path.dirname(os.path.abspath(stru_path))
            valences = []
            try:
                for elem in species_order:
                    if elem not in pp_map or pp_map[elem] in (None, ""):
                        raise KeyError(f"pp 映射缺少元素 {elem} 的赝势文件名")
                    pp_name_or_path = pp_map[elem]
                    upf_path = pp_name_or_path if os.path.isabs(pp_name_or_path) else os.path.normpath(os.path.join(base_dir, pp_name_or_path))
                    if not os.path.isfile(upf_path):
                        raise FileNotFoundError(f"赝势文件不存在：{upf_path}")
                    z = _parse_upf_valence(upf_path)
                    valences.append(int(round(z)))
                valence_line = ' '.join(map(str, valences))
            except Exception as e:
                warnings.warn(f"从 pp 解析 valence 失败：{e}。将使用占位提示继续运行。")
                use_placeholder = True
                valence_line = placeholder_valence_line


    # —— 仿写你的 generate_input_optical 的拼接风格 ——
    input_text += f"\nPOLARIZATION\n"
    input_text += f"{{\n"
    input_text += f"    occ_band       {n_occu}\n"
    input_text += f"    nk1            {nk1}\n"
    input_text += f"    nk2            {nk2}\n"
    input_text += f"    nk3            {nk3}\n"
    input_text += f"    atom_type      {atom_type}\n"
    input_text += f"    stru_file      {stru_file_field}\n"
    input_text += f"    valence_e      {valence_line}\n"
    input_text += f"}}\n"

    return input_text

def generate_input_anc(input_text, dim, lattice_vectors, method, mp_density, energy_range):

    e_fermi = 0
    # 检查energy_range的长度，来决定是直接使用还是基于费米能级计算

    e_low, e_high = e_fermi + energy_range[0], e_fermi + energy_range[1]


    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, mp_density)  # AHC 默认密度 0.008

    input_text += f"\nANC\n"
    input_text += f"{{\n"
    input_text += f"    method                    {method}\n"
    input_text += f"    fermi_range               {e_low} {e_high}\n"
    input_text += f"    de                        0.01\n"
    input_text += f"    eta                       0.10\n"
    input_text += f"    integrate_grid            {' '.join(map(str, integrate_grid))}\n"
    input_text += f"}}\n"
    # print(input_text)
    return input_text



def generate_input_wilsonloop(input_text, n_occu, occu_switch, dim, lattice_vectors, method):
    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.05)  # Wilson Loop 密度

    input_text += f"\nWILSON_LOOP\n"
    input_text += f"{{\n"
    if occu_switch == 0:
        input_text += f"    occ_band                  {n_occu}\n"  # autoset to VBM
    elif occu_switch > 0:
        input_text += f"    occ_band                  {occu_switch}\n"  # set to specific band
    else:
        pass
    input_text += f"    k_start  0.0  0.0  0.5\n"
    input_text += f"    k_vect1  1.0  0.0  0.0\n"
    input_text += f"    k_vect2  0.0  0.5  0.0\n"
    input_text += f"    nk1  101\n"
    input_text += f"    nk2  101\n"
    input_text += f"}}\n"
    return input_text


def generate_input_chern(input_text, n_occu, occu_switch, dim, lattice_vectors, method, mp_density):
    integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, mp_density)  # Chern 密度

    input_text += f"\nCHERN_NUMBER\n"
    input_text += f"{{\n"
    input_text += f"    method                    {method}\n"
    if occu_switch == 0:
        input_text += f"    occ_band                  {n_occu}\n"  # autoset to VBM
    elif occu_switch > 0:
        input_text += f"    occ_band                  {occu_switch}\n"  # set to specific band
    else:
        pass
    input_text += f"    integrate_mode            Grid\n"
    input_text += f"    integrate_grid            {' '.join(map(str, integrate_grid))}\n"
    # input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
    # input_text += f"    adaptive_grid_threshold   1000\n"
    input_text += f"    k_start  0 0 0\n"
    input_text += f"    k_vect1  1 0 0\n"
    input_text += f"    k_vect2  0 1 0\n"
    input_text += f"}}\n"
    return input_text


def generate_input_berry(input_text, ase_stru: Atoms, n_occu, occu_switch, kline_density, dim, kmode, tolerance, knum, kpath, method, mp_density):
    lattice_vectors = ase_stru.get_cell()

    input_text += f"\nBERRY_CURVATURE\n"
    input_text += f"{{\n"
    input_text += f"    method                    {method}\n"
    if occu_switch == 0:
        input_text += f"    occ_band                  {n_occu}\n"  # autoset to VBM
    elif occu_switch > 0:
        input_text += f"    occ_band                  {occu_switch}\n"  # set to specific band
    else:
        pass

    if kmode in ["mp", "mesh"]:
        integrate_grid, adaptive_grid, pbc = get_k_mesh_from_dim(lattice_vectors, dim, 0.02)  # berry mesh 默认密度
        input_text += f"    kpoint_mode               mp\n"
        input_text += f"    mp_grid                   {' '.join(map(str, integrate_grid))}\n"
        input_text += f"    adaptive_grid             {' '.join(map(str, adaptive_grid))}\n"
        input_text += f"    adaptive_grid_threshold   1000\n"
    else:
        band_line,  kpoint_label_line, kpoint_num_in_line = kpath_generator(ase_stru, kline_density, dim, tolerance=tolerance, knum=knum, kpath=kpath)
        num_lines = len(band_line) - 1
        input_text += f"    kpoint_mode               {kmode}\n"
        input_text += f"    kpoint_num                {num_lines}\n"
        input_text += f"    high_symmetry_kpoint\n"
        for line in band_line:
            input_text += f"   {line}\n"
    input_text += f"}}\n"
    return input_text


def extract_data_from_log(out_suffix_path):
    """
    从 running_scf.log 文件中提取数据（增强版）。

    规则：
    - 常规情况下：从全文件解析 e_tot、e_fermi、n_occu、n_bands、n_elec（遇到多次取最后一次）。
    - 若出现 'nelec_delta is NOT zero' 行：从该行 **之后** 的内容重新解析并覆盖 n_elec / n_occu / n_bands。
      * n_elec 优先取 'nelec now'
      * 若没有 'nelec now'，再退而求其次找 'AUTOSET number of electrons'
      * 不会把 'total electron number of element X' 当作总电子数

    返回：
        e_tot (float): 总能量
        e_fermi (float): 费米能级
        n_occu (int): 占据态数目（occupied bands）
        n_bands (int): 总带数目（NBANDS）
        n_elec (int): 电子数目（首选 'nelec now'）
    """
    log_file_path = os.path.join(out_suffix_path, "running_scf.log")
    if not os.path.isfile(log_file_path):
        print(f"错误：未找到日志文件 {log_file_path}")
        return

    # 结果变量（常规/默认解析结果）
    e_tot = None
    e_fermi = None
    n_occu = None
    n_bands = None
    n_elec = None

    # “nelec_delta 警告之后”的覆盖结果（若触发，将用这些覆盖默认结果）
    post_occu = None
    post_bands = None
    post_elec = None

    # 标记：是否进入“nelec_delta 警告之后”的区域
    after_nelec_delta = False

    # 预编译正则
    re_etot   = re.compile(r"!FINAL_ETOT_IS\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)")
    re_fermi  = re.compile(r"EFERMI\s*=\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)")
    re_occu   = re.compile(r"occupied\s+bands\s*=\s*(\d+)", re.IGNORECASE)
    re_bands  = re.compile(r"NBANDS\s*=\s*(\d+)", re.IGNORECASE)

    # 不同来源的电子数（注意优先级：nelec now > AUT0SET ... > number of electrons）
    re_nelec_now   = re.compile(r"nelec\s+now\s*:\s*=\s*(\d+)", re.IGNORECASE)
    re_autoset_ele = re.compile(r"AUTOSET\s+number\s+of\s+electrons\s*:\s*=\s*(\d+)", re.IGNORECASE)
    re_num_elec    = re.compile(r"number\s+of\s+electrons\s*=?\s*(:|=)?\s*([+-]?\d+(?:\.\d+)?)", re.IGNORECASE)

    # 警告行
    re_nelec_delta_warn = re.compile(r"nelec_delta\s+is\s+NOT\s+zero", re.IGNORECASE)

    # 非总电子数的提示（仅记录，不计入 n_elec）
    re_element_elec = re.compile(r"total\s+electron\s+number\s+of\s+element\s+\S+\s*=\s*\d+", re.IGNORECASE)

    try:
        with open(log_file_path, "r", encoding="utf-8", errors="ignore") as log_file:
            for raw_line in log_file:
                line = raw_line.strip()

                # 全局解析（e_tot/e_fermi 始终取最后一次）
                m = re_etot.search(line)
                if m:
                    e_tot = float(m.group(1))
                    continue

                m = re_fermi.search(line)
                if m:
                    e_fermi = float(m.group(1))
                    continue

                # 进入“警告之后”模式
                if re_nelec_delta_warn.search(line):
                    after_nelec_delta = True
                    # 进入该模式后，将使用 post_* 收集覆盖项
                    # 不在此处返回，继续向下扫描“之后”的值
                    continue

                # occupied bands / NBANDS / n_elec 的解析
                # 根据是否在“警告之后”写入不同变量
                if after_nelec_delta:
                    m = re_occu.search(line)
                    if m:
                        post_occu = int(m.group(1))
                        continue

                    m = re_bands.search(line)
                    if m:
                        post_bands = int(m.group(1))
                        continue

                    # 电子数优先级：nelec now > AUTOSET ...
                    m = re_nelec_now.search(line)
                    if m:
                        post_elec = int(m.group(1))
                        continue
                    m = re_autoset_ele.search(line)
                    if m:
                        # 只有在 post_elec 尚未由 'nelec now' 设定时才用它
                        if post_elec is None:
                            post_elec = int(m.group(1))
                        continue
                    m = re_num_elec.search(line)
                    if m:
                        # 尽量避免把“每元素总电子数”等误判为总电子数
                        if not re_element_elec.search(line):
                            # 只有当既没有 nelec now 也没有 AUTOSET 时，才兜底用它
                            if post_elec is None:
                                try:
                                    post_elec = int(float(m.group(2)))
                                except ValueError:
                                    pass
                        continue
                else:
                    # 常规模式（在警告之前）
                    m = re_occu.search(line)
                    if m:
                        n_occu = int(m.group(1))
                        continue

                    m = re_bands.search(line)
                    if m:
                        n_bands = int(m.group(1))
                        continue

                    # 电子数：先记录，但如果后续有“警告之后”的值会被覆盖
                    m = re_nelec_now.search(line)
                    if m:
                        n_elec = int(m.group(1))
                        continue
                    m = re_autoset_ele.search(line)
                    if m:
                        # 若未由 nelec now 设定，则可先用
                        if n_elec is None:
                            n_elec = int(m.group(1))
                        continue
                    m = re_num_elec.search(line)
                    if m:
                        if not re_element_elec.search(line):
                            try:
                                n_elec = int(float(m.group(2)))
                            except ValueError:
                                pass
                        continue

        # 如果出现了警告，则使用“警告之后”的覆盖值
        if after_nelec_delta:
            if post_occu is not None:
                n_occu = post_occu
            if post_bands is not None:
                n_bands = post_bands
            if post_elec is not None:
                n_elec = post_elec

        # 基本检查
        if e_tot is None or e_fermi is None:
            print("错误：无法从 running_scf.log 中提取 e_tot 或 e_fermi")
            return

        print(f"{log_file_path} 提取完成。")
        print(f"E_TOTAL (eV)   =  {e_tot}")
        print(f"E_FERMI (eV)   =  {e_fermi}")
        print(f"OCCUPIED_BANDS =  {n_occu}")
        print(f"NBANDS         =  {n_bands}")
        print(f"N_ELEC         =  {n_elec}")

        return e_tot, e_fermi, n_occu, n_bands, n_elec

    except Exception as e:
        print(f"读取/解析日志时发生异常：{e}")
        return


def stru_pp_orb(pseudo_dir, orbital_dir, f_stru, ase_stru: Atoms, f_out=None):
    """
    Update the pseudopotential and orbital filenames in the structure file.

    Args:
        pseudo_dir (str): Directory path for the pseudopotential files.
        orbital_dir (str): Directory path for the orbital files.
        f_stru (str): Path to the input structure file.
        ase_stru (ase.Atoms): ASE Atoms object representing the structure.
        f_out (str, optional): Path to the output structure file. If not provided, the input file will be overwritten.

    Returns:
        None
    """
    with open(f_stru, 'r') as file:
        content = file.read()

    atoms_pp = ase_stru.info['pp']
    atoms_orb = ase_stru.info['basis']
    atoms_all = ase_stru.get_chemical_symbols()
    atoms_dict = {}
    for idx, atoms_all_name in enumerate(atoms_all):
        if atoms_all_name not in atoms_dict:
            atoms_dict[atoms_all_name] = [ase_stru.get_masses()[idx], atoms_pp[atoms_all_name], atoms_orb[atoms_all_name]]

    for key, values in atoms_dict.items():
        element = key
        mass = values[0]
        pp = values[1]
        orb = values[2]
        pp_filename = os.path.join(pseudo_dir, pp)
        orb_filename = os.path.join(orbital_dir, orb)
        content = re.sub(rf"{element}\s+\d+\.\d+\s+.*\.upf", f"{element} {mass} {pp_filename}", content)
        content = re.sub(f".*{element}_.+\\.orb", f"{orb_filename}", content)

    with open(f_out, 'w') as file:
        file.write(content)
    return




def main():
    # 创建ArgumentParser对象
    parser = argparse.ArgumentParser(description='Input Generator Script.')
    parser.add_argument('-i', '--input', '--in', type=str, default='./', help='Set ABACUS SCF running directory, DEFAULT is ./, i.e. PWD')
    parser.add_argument('-o', '--output', '--out', type=str, default=None, help='Set output directory to run PYATB. DEFAULT is to make a new directory in the input_scf_dir/pyatb')
    parser.add_argument('--band', action='store_true', help='Band Structure calculation')
    parser.add_argument('--kline', type=float, default=0.01, help='Density of K-path, defalut is 0.01 2*pi/Angstrom.')
    parser.add_argument('--knum', type=int, default=0, help='Number oof kpoints on each line of K-path, defalut is 0.')
    parser.add_argument('--kpath', type=str, default=None, help='Line of K-path, you can use GKMG to perform a 2D K-path, defalut is None')
    parser.add_argument('--kmode', type=str, default='line', help='Mode of K-path, defalut is line, you can choose mp or line.')
    parser.add_argument('--dim', type=str, default='3', help='Band K-path dimension. Default value is 3 which represents the 3D band path, when using the "1, 1, 0" format of the string represents the generation of only ab direction k point')
    parser.add_argument('--tolerance', '--tol', type=float, default=1e-3, help='Tolerance for determining Bravais lattice. Default is 1e-3.')
    parser.add_argument('--pdos', action='store_true', help='Projected DOS calculation, defalut energy range referred to Fermi level is 4 which means [-4, 4], you can also set as --erange="-6 8" to set [-6, 8]')
    parser.add_argument('--findnodes', '--fnodes', action='store_true', help='Find Weyl nodes calculation')
    parser.add_argument('--erange', type=str, default='8', help='Energy range referred to Fermi level, defalut is 8 means [-8, 8], you can also set as --erange="-6 8" to set [-6, 8]')
    parser.add_argument('--frange', type=float, default=1, help='Fermi energy range referred to the Fermi level, defalut is -1 1')
    parser.add_argument('--optical', action='store_true', help='Optical conductivity with automatically read band occupation numbers')
    parser.add_argument('--jdos', action='store_true', help='JDOS with automatically read band occupation numbers')
    parser.add_argument('--shift','--shift_current', action='store_true', help='Shift current with automatically read band occupation numbers')
    parser.add_argument('--polar','--polarization', action='store_true', help='Calculate the polarization of periodic solids by Modern Theory of Polarization, namely Berry phase theory')
    parser.add_argument('--valence', '--valence_e', type=str, default='auto', help='for polarization calculation, default is auto, you can set like "14 3 4" for each of 3 elements.')
    parser.add_argument('--orange', type=float, nargs="+", default=[0, 10], help='Optical frequency range, default is [0, 10] eV')
    parser.add_argument('--mp','--density', type=float, default=0.10, help='MP Grid density along every reciprocal vector, default is 0.10 2*pi/Angstrom.')
    parser.add_argument('--bandunfolding','--unfolding','--bandunfold','--unfold', action='store_true', help='Band Structure Unfolding calculation, use --bandrange to set the band range, and --matrix to set the unfolding matrix.')
    parser.add_argument('--m_matrix', '--matrix', type=str, default='1 0 0  0 1 0  0 0 1', help='Band unfolding matrix in the Input file.')
    parser.add_argument('--ahc', action='store_true', help='Anomalous Hall Conductivity with defalut setting Input generated.')
    parser.add_argument('--anc', action='store_true', help='Anomalous Nerst Conductivity, you can modify the Input file by --erange.')
    parser.add_argument('--berry', '-b',action='store_true', help='Berry Curvature with defalut setting Input generated. Please modify the Input file manually.')
    parser.add_argument('--occu', type=int, default=0, help='Berry Curvature occupied band switch. -1 means no occupied band, 0 means autoset the VBM, 66 number means the 66th band, etc.')
    parser.add_argument('--method', type=int, default=0, help='Method for calculating Berry curvature. 0 means direct calculation, 1 means calculation by Kubo formula.')
    parser.add_argument('--bcd','--berry_curvature_dipole', action='store_true', help='BERRY_CURVATURE_DIPOLE with defalut setting Input generated.')
    parser.add_argument('--cpge', action='store_true', help='CPGE with defalut setting Input generated. Please modify the Input file manually.')
    parser.add_argument('--chern', action='store_true', help='Chern Number with defalut setting Input generated.')
    parser.add_argument('--wilson','--wl','--wilson_loop', action='store_true', help='WILSON_LOOP with defalut setting Input generated.')
    parser.add_argument('--fatband','--pband','--projectedband', action='store_true', help='Projected Band Structures with defalut setting Input generated.')
    parser.add_argument('--spintexture','--spin','--spintex', action='store_true', help='Spin Texture Input generated. DEFAULT of Kpoint mode is line, you can set --kmode=mp to use MP grid.')
    parser.add_argument('--bandrange','--brange', type=str, default=None, help='Band range for FATBAND or BAND UNFOLDING, defalut is +- 50 bands about the occupied band.')
    parser.add_argument('--max_kpoint_num','--maxkpt',type=int, default=2000, help='Max parallel kpoint number in one iteration.')
    parser.add_argument('--fs','--fermisurface', action='store_true', help='Fermi surface with defalut setting Input generated.')
    calculator = parser.add_mutually_exclusive_group()
    calculator.add_argument('--abacus', action='store_true', help='Use ABACUS for SCF calculator.')
    calculator.add_argument('--hamgnn', action='store_true', help='Use HamGNN for .csr generation.')

    args = parser.parse_args()

    # 默认是 abacus，如果都没指定，就手动设置 abacus=True
    if not args.abacus and not args.hamgnn:
        args.abacus = True

    # 现在你可以根据选择来决定计算器
    if args.abacus:
        print("Using ABACUS as calculator.")
    elif args.hamgnn:
        print("Using HamGNN as calculator.")




    if args.hamgnn:
        directory_path =  args.input
        directory_path = os.path.abspath(directory_path)

        output_pyatb_path =  args.output

        out_suffix_path = directory_path

        properties_json_file = os.path.join(out_suffix_path, "properties.json")
        latname = None
        nspin = 1 # 先写死，后面替换
        # 使用get方法从字典中获取键值，如果键不存在，则返回空字符串 ''
        pp_dir  = None
        orb_dir = None

        # e_tot, e_fermi, n_occu, n_bands, n_elec = extract_data_from_log(out_suffix_path)

        max_val, e_fermi, band_gap, n_elec, n_occu = extract_properties_from_json(properties_json_file)

        n_bands = int(n_elec) # 先这样设置

        stru_file_path = os.path.join(directory_path, "STRU")
        with open(stru_file_path, 'r') as f_s:
            if latname == 'none':
                i_latname = None
            else:
                i_latname = latname
            ase_stru = read_abacus_stru(f_s, i_latname, True)
        lattice_constant = 1.0 # unit is Angstrom
        lattice_vectors = ase_stru.get_cell()


        input_text = generate_hamgnn_input_init(nspin, e_fermi, out_suffix_path, lattice_constant, lattice_vectors, args.max_kpoint_num)
        
        current_directory = os.path.abspath(os.getcwd())

        # 判断是否提供了output_pyatb_path，如果提供了，那么就在目录下生成相应名称的文件夹
        if output_pyatb_path is not None:
            pyatb_directory = os.path.join(directory_path, output_pyatb_path)
            os.makedirs(pyatb_directory, exist_ok=True)
            shutil.copy(os.path.join(current_directory, "STRU"), os.path.join(pyatb_directory, "STRU"))
        # 判断当前目录是否在 directory_path 的子目录中，且不是 directory_path 本身
        elif os.path.commonpath([current_directory, directory_path]) == directory_path and current_directory != directory_path:
            # 如果当前目录是在 directory_path 内部（但不等于 directory_path 本身）
            pyatb_directory = current_directory
        else:
            # 如果当前目录等于 directory_path 或者不在 directory_path 内部
            pyatb_directory = os.path.join(directory_path, "pyatb")
            os.makedirs(pyatb_directory, exist_ok=True)

        get_e_file_path = os.path.join(pyatb_directory, "get_Energy_HamGNN.out")
        with open(get_e_file_path, "w") as output_file:
            output_file.write(f"E_TOTAL (eV)   =  None_yet\n")
            output_file.write(f"E_FERMI (eV)   =  {e_fermi}\n")
            output_file.write(f"BAND_GAP       =  {band_gap}\n")
            output_file.write(f"NBANDS         =  None_yet \n")
            output_file.write(f"NELEC          =  {n_elec}\n")
            output_file.write(f"Occupied bands =  {n_occu} # estimated by NELEC/2 \n")

        if args.bandrange:
            if ' ' in args.bandrange:
                band_range = [int(band) for band in args.bandrange.split()]
            elif len(args.bandrange.split()) == 1:
                band_range = int(args.bandrange.split()[0])
                band_range = [max(1, n_occu - band_range), min(n_bands, n_occu + band_range)]  # 上下band_range条
        else:   
            band_range = [max(1, n_occu - 100), min(n_bands, n_occu + 100)]  # 上下100条， 如果最小值小于1那就等于1，如果最大值大于n_bands，那就等于n_bands

        if args.orange:
            if len(args.orange) == 1:
                omega_range = [0.0, args.orange[0]]
            elif len(args.orange) == 2:
                omega_range = [args.orange[0], args.orange[1]]
        
        if args.erange:
            erange_value = args.erange.strip()  # 去除前后空格
            # 判断是单一值还是范围
            if ' ' in erange_value:
                # 如果是两个数值的字符串，按空格拆分并转换为浮动数值列表
                energy_range = [float(v) for v in erange_value.split()]
            else:
                # 如果是单一值，转换为[-value, value]的范围
                num = float(erange_value)
                energy_range = [-num, num]

        # 传递目录路径到每个函数
        if args.band:
            input_text = generate_input_band(input_text, ase_stru, args.kline, args.dim,  args.kmode, args.tolerance,  args.knum,  args.kpath)
        if args.fatband:
            input_text = generate_input_fatband(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath,  band_range, n_bands)
        if args.spintexture:
            input_text = generate_input_spintexture(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath,  band_range, n_bands)
        if args.bandunfolding:
            input_text = generate_input_bandunfold(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.m_matrix, args.kpath,  band_range, n_bands)
        if args.ahc:
            input_text = generate_input_ahc(input_text,  args.dim, lattice_vectors, args.mp)
        if args.anc:
            input_text = generate_input_anc(input_text,  args.dim, lattice_vectors, args.method, args.mp, energy_range)
        if args.berry:
            input_text = generate_input_berry(input_text, ase_stru, n_occu, args.occu, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath, args.method, args.mp)
        if args.bcd:
            input_text = generate_input_bcd(input_text,  args.dim, lattice_vectors,  e_fermi, energy_range)
        if args.cpge:
            input_text = generate_input_bcd(input_text,  args.dim,  lattice_vectors, e_fermi,energy_range)
        if args.pdos:
            input_text = generate_input_pdos(input_text, args.dim, lattice_vectors, e_fermi, energy_range)
        if args.findnodes:
            input_text = generate_input_findnodes(input_text, args.dim, lattice_vectors, e_fermi, energy_range)
        if args.optical:
            input_text = generate_input_optical(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
        if args.polar:
            input_text = generate_input_polarization(input_text, args.dim, lattice_vectors, n_occu, ase_stru, stru_file_path, args.valence, args.mp, stru_file_field="STRU")
        if args.jdos:
            input_text = generate_input_jdos(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
        if args.shift:
            input_text = generate_input_shift(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
        if args.chern:
            input_text = generate_input_chern(input_text,  n_occu, args.occu, args.dim, lattice_vectors, args.method, args.mp)
        if args.wilson:
            input_text = generate_input_wilsonloop(input_text,  n_occu, args.occu, args.dim, lattice_vectors, args.method)


        input_file_path = os.path.join(pyatb_directory, "Input")
        with open(input_file_path, "w") as file:
            file.write(input_text)

        print(f"[PYATB with HamGNN] Input and STRU files for PyATB is generated in: {pyatb_directory}")

        return

    



    if args.abacus:
        directory_path =  args.input
        output_pyatb_path =  args.output

        if os.path.isdir(directory_path):
            input_file_path = os.path.join(directory_path, "INPUT")
        elif os.path.isfile(directory_path):
            input_file_path = directory_path
            directory_path = os.path.dirname(directory_path)
        else:
            print(f"[PYATB with ABACUS] 路径或文件 {directory_path} 不存在。")
            return

        # 把directory_path转换为绝对路径
        directory_path = os.path.abspath(directory_path)

        variables_dict_init = parse_input_file(input_file_path)
        if variables_dict_init:
            if "suffix" in variables_dict_init:
                out_suffix_path = os.path.join(directory_path, "OUT." + variables_dict_init["suffix"])
            else:
                out_suffix_path = os.path.join(directory_path, "OUT.ABACUS")

            input_file_full = os.path.join(out_suffix_path, "INPUT")
            variables_dict_full = parse_input_file(input_file_full)
            latname = variables_dict_full["latname"]
            nspin = int(variables_dict_full["nspin"])
            # 使用get方法从字典中获取键值，如果键不存在，则返回空字符串 ''
            pp_dir = os.path.join(directory_path, variables_dict_full.get("pseudo_dir", ''))
            orb_dir = os.path.join(directory_path, variables_dict_full.get("orbital_dir", ''))

            e_tot, e_fermi, n_occu, n_bands, n_elec = extract_data_from_log(out_suffix_path)

            stru_file_path = os.path.join(directory_path, "STRU")
            with open(stru_file_path, 'r') as f_s:
                if latname == 'none':
                    i_latname = None
                else:
                    i_latname = latname
                ase_stru = read_abacus_stru(f_s, i_latname, True)
            lattice_constant = 1.0 # unit is Angstrom
            lattice_vectors = ase_stru.get_cell()

            input_text = generate_input_init(nspin, e_fermi, out_suffix_path, lattice_constant, lattice_vectors, args.max_kpoint_num)
            
            current_directory = os.path.abspath(os.getcwd())
            # 判断是否提供了output_pyatb_path，如果提供了，那么就在目录下生成相应名称的文件夹
            if output_pyatb_path is not None:
                pyatb_directory = os.path.join(directory_path, output_pyatb_path)
                os.makedirs(pyatb_directory, exist_ok=True)
            # 判断当前目录是否在 directory_path 的子目录中，且不是 directory_path 本身
            elif os.path.commonpath([current_directory, directory_path]) == directory_path and current_directory != directory_path:
                # 如果当前目录是在 directory_path 内部（但不等于 directory_path 本身）
                pyatb_directory = current_directory
            else:
                # 如果当前目录等于 directory_path 或者不在 directory_path 内部
                pyatb_directory = os.path.join(directory_path, "pyatb")
                os.makedirs(pyatb_directory, exist_ok=True)
            get_e_file_path = os.path.join(pyatb_directory, "get_Energy.out")
            with open(get_e_file_path, "w") as output_file:
                output_file.write(f"E_TOTAL (eV)   =  {e_tot}\n")
                output_file.write(f"E_FERMI (eV)   =  {e_fermi}\n")
                output_file.write(f"Occupied bands =  {n_occu}\n")
                output_file.write(f"NBANDS         =  {n_bands}\n")
                output_file.write(f"NELEC          =  {n_elec}\n")

            # 把 stru_file_path 生成一份到pyatb文件夹下，注意文件夹路径改变，PP和ORB的位置可能改变，此处进行判断和修改
            if pp_dir is None and orb_dir is None: # 如果没有找到pseudo_dir和orbital_dir，不修改STRU文件
                pass
            else:
                stru_pp_orb(pp_dir, orb_dir, f_stru=stru_file_path, ase_stru=ase_stru, f_out=os.path.join(pyatb_directory, "STRU"))


            if args.bandrange:
                if ' ' in args.bandrange:
                    band_range = [int(band) for band in args.bandrange.split()]
                elif len(args.bandrange.split()) == 1:
                    band_range = int(args.bandrange.split()[0])
                    band_range = [max(1, n_occu - band_range), min(n_bands, n_occu + band_range)]  # 上下band_range条
            else:   
                band_range = [max(1, n_occu - 100), min(n_bands, n_occu + 100)]  # 上下100条， 如果最小值小于1那就等于1，如果最大值大于n_bands，那就等于n_bands

            if args.orange:
                if len(args.orange) == 1:
                    omega_range = [0.0, args.orange[0]]
                elif len(args.orange) == 2:
                    omega_range = [args.orange[0], args.orange[1]]
            
            if args.erange:
                erange_value = args.erange.strip()  # 去除前后空格
                # 判断是单一值还是范围
                if ' ' in erange_value:
                    # 如果是两个数值的字符串，按空格拆分并转换为浮动数值列表
                    energy_range = [float(v) for v in erange_value.split()]
                else:
                    # 如果是单一值，转换为[-value, value]的范围
                    num = float(erange_value)
                    energy_range = [-num, num]

            # 传递目录路径到每个函数
            if args.band:
                input_text = generate_input_band(input_text, ase_stru, args.kline, args.dim,  args.kmode, args.tolerance,  args.knum,  args.kpath)
            if args.fatband:
                input_text = generate_input_fatband(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath,  band_range, n_bands)
            if args.spintexture:
                input_text = generate_input_spintexture(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath,  band_range, n_bands)
            if args.bandunfolding:
                input_text = generate_input_bandunfold(input_text, ase_stru, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.m_matrix, args.kpath,  band_range, n_bands)
            if args.ahc:
                input_text = generate_input_ahc(input_text,  args.dim, lattice_vectors, args.mp)
            if args.anc:
                input_text = generate_input_anc(input_text,  args.dim, lattice_vectors, args.method, args.mp, energy_range)
            if args.berry:
                input_text = generate_input_berry(input_text, ase_stru, n_occu, args.occu, args.kline, args.dim, args.kmode, args.tolerance, args.knum, args.kpath, args.method, args.mp)
            if args.bcd:
                input_text = generate_input_bcd(input_text,  args.dim, lattice_vectors,  e_fermi, energy_range)
            if args.cpge:
                input_text = generate_input_bcd(input_text,  args.dim,  lattice_vectors, e_fermi,energy_range)
            if args.pdos:
                input_text = generate_input_pdos(input_text, args.dim, lattice_vectors, e_fermi, energy_range)
            if args.findnodes:
                input_text = generate_input_findnodes(input_text, args.dim, lattice_vectors, e_fermi, energy_range)
            if args.optical:
                input_text = generate_input_optical(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
            if args.polar:
                input_text = generate_input_polarization(input_text, args.dim, lattice_vectors, n_occu, ase_stru, stru_file_path, args.valence, args.mp, stru_file_field="STRU")
            if args.jdos:
                input_text = generate_input_jdos(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
            if args.shift:
                input_text = generate_input_shift(input_text, args.dim, lattice_vectors, n_occu, omega_range, args.mp)
            if args.chern:
                input_text = generate_input_chern(input_text,  n_occu, args.occu, args.dim, lattice_vectors, args.method, args.mp)
            if args.wilson:
                input_text = generate_input_wilsonloop(input_text,  n_occu, args.occu, args.dim, lattice_vectors, args.method)


            input_file_path = os.path.join(pyatb_directory, "Input")
            with open(input_file_path, "w") as file:
                file.write(input_text)

            print(f"[PYATB with ABACUS]Input and STRU files for PyATB is generated in: {pyatb_directory}")
            return
        else:
            print("[PYATB with ABACUS] 未能从INPUT文件中解析出任何变量。")
    

if __name__ == "__main__":
    main()
