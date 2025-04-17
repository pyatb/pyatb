from pyatb import RANK, INPUT_PATH, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.fermi.fat_band import Fat_Band
from pyatb.io import read_input
from pyatb.tb import tb
from pyatb.tb import multiXR
from pyatb.io.abacus_read_xr import abacus_readHR, abacus_readSR, abacus_readrR
from pyatb.io.wannier90_read_tb import wannier90_readHR, wannier90_readTB
from pyatb.fermi import *  
from pyatb.berry import *
from pyatb.transport import *
import os

def main():
    # get user INPUT
    INPUT, function_switch, bool_need_rR = read_input(os.path.join(INPUT_PATH, 'Input'))

    input_parameters = INPUT['INPUT_PARAMETERS']
    nspin = input_parameters['nspin']
    lattice_constant = INPUT['LATTICE']['lattice_constant']
    lattice_vector = INPUT['LATTICE']['lattice_vector']
    package = input_parameters['package']
    sparse_format = input_parameters['sparse_format']
    max_kpoint_num = input_parameters['max_kpoint_num']

    # initialize the tb class
    m_tb = tb(nspin, lattice_constant, lattice_vector, max_kpoint_num)

    if package == 'ABACUS':
        if nspin != 2:
            HR = abacus_readHR(**input_parameters)
        else:
            HR_up_route = input_parameters['HR_route'][0]
            HR_dn_route = input_parameters['HR_route'][1]
            HR_unit = input_parameters['HR_unit']

            HR_up = abacus_readHR(nspin, HR_up_route, HR_unit)
            HR_dn = abacus_readHR(nspin, HR_dn_route, HR_unit)

        SR = abacus_readSR(**input_parameters)

        if bool_need_rR:
            rR = abacus_readrR(**input_parameters)
    elif package == 'WANNIER90':
        if nspin == 2:
            raise ValueError('WANNIER90 only for nspin = 1 or 4 !')

        if input_parameters['w90_TB_has_r']:
            HR, SR, *rR = wannier90_readTB(**input_parameters)
        else:
            HR, SR = wannier90_readHR(**input_parameters)
    
    # set tb
    if sparse_format == 0:
        isSparse = False
    else:
        isSparse = True

    if nspin != 2:
        m_tb.set_solver_HSR(HR, SR, isSparse)
    else:
        m_tb.set_solver_HSR_spin2(HR_up, HR_dn, SR, isSparse)

    if bool_need_rR:
        m_tb.set_solver_rR(rR[0], rR[1], rR[2], isSparse)

    # output information
    if RANK == 0:
        with open(RUNNING_LOG, 'a') as f:
            f.write('\n')
            f.write('\n------------------------------------------------------')
            f.write('\n|                                                    |')
            f.write('\n|                     Structure                      |')
            f.write('\n|                                                    |')
            f.write('\n------------------------------------------------------')
            f.write('\n\n')

            f.write('Lattice Constant : %.6f Angstrom\n' %(lattice_constant))
            f.write('Lattice Vector   : \n')
            f.write('%15.6f %15.6f %15.6f \n%15.6f %15.6f %15.6f \n%15.6f %15.6f %15.6f \n'%(
                lattice_vector[0, 0], lattice_vector[0, 1], lattice_vector[0, 2], 
                lattice_vector[1, 0], lattice_vector[1, 1], lattice_vector[1, 2],
                lattice_vector[2, 0], lattice_vector[2, 1], lattice_vector[2, 2])
            )

    # calculation
    if function_switch['FERMI_ENERGY']:
        fermi_energy_parameters = INPUT['FERMI_ENERGY']
        cal_FE = Fermi_Energy(m_tb)
        input_parameters['fermi_energy'] = cal_FE.calculate_fermi_energy(**fermi_energy_parameters, **input_parameters)

    if function_switch['BAND_STRUCTURE']:
        band_structure_parameters = INPUT['BAND_STRUCTURE']
        wf_collect = band_structure_parameters['wf_collect']
        cal_BS = Band_Structure(m_tb, wf_collect)
        fermi_energy = input_parameters['fermi_energy']
        cal_BS.calculate_band_structure(fermi_energy=fermi_energy, **band_structure_parameters)
        if RANK == 0:
            cal_BS.print_plot_script()

    if function_switch['BANDUNFOLDING']:
        bandunfolding_parameters = INPUT['BANDUNFOLDING']
        M_matrix = bandunfolding_parameters['m_matrix'].reshape(3, 3)
        bandunfolding_parameters['m_matrix'] = M_matrix
        cal_BUNF = Bandunfolding(m_tb)
        cal_BUNF.calculate_bandunfolding(**bandunfolding_parameters)
        if RANK == 0:
            fermi_energy = input_parameters['fermi_energy']
            cal_BUNF.print_plot_script(fermi_energy)

    if function_switch['BANDUNFOLDING_SPIN_TEXTURE']:
        bandunfolding_spin_texture_parameters = INPUT['BANDUNFOLDING_SPIN_TEXTURE']
        M_matrix = bandunfolding_spin_texture_parameters['m_matrix'].reshape(3, 3)
        bandunfolding_spin_texture_parameters['m_matrix'] = M_matrix
        cal_BUNF_ST = Bandunfolding_Spin_Texture(m_tb)
        cal_BUNF_ST.calculate_bandunfolding_spin_texture(**bandunfolding_spin_texture_parameters)
        if RANK == 0:
            fermi_energy = input_parameters['fermi_energy']
            cal_BUNF_ST.print_plot_script(fermi_energy)

    if function_switch['FAT_BAND']:
        fatband_parameters = INPUT['FAT_BAND']
        cal_FAT = Fat_Band(m_tb)
        cal_FAT.calculate_fatband(**fatband_parameters)
        if RANK == 0:
            fermi_energy = input_parameters['fermi_energy']
            cal_FAT.print_plot_script(fermi_energy)

    if function_switch['FERMI_SURFACE']:
        fermi_surface_parameters = INPUT['FERMI_SURFACE']
        cal_FS = Fermi_Surface(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_FS.calculate_fermi_surface(fermi_energy=fermi_energy, **fermi_surface_parameters)
        if RANK == 0:
            cal_FS.print_plot_script()

    if function_switch['FIND_NODES']:
        find_nodes_parameters = INPUT['FIND_NODES']
        cal_nodes = Find_Nodes(m_tb)
        cal_nodes.calculate_nodes(**find_nodes_parameters)
        if RANK == 0:
            cal_nodes.print_plot_script()

    if function_switch['JDOS']:
        jdos_parameters = INPUT['JDOS']
        cal_JDOS = JDOS(m_tb)
        cal_JDOS.calculate_jdos(**jdos_parameters)
        if RANK == 0:
            cal_JDOS.print_plot_script()

    if function_switch['PDOS']:
        pdos_parameters = INPUT['PDOS']
        cal_PDOS = PDOS(m_tb)
        cal_PDOS.calculate_dos(**pdos_parameters)
        if RANK == 0:
            fermi_energy = input_parameters['fermi_energy']
            e_range = pdos_parameters['e_range']
            cal_PDOS.print_plot_script(fermi_energy, e_range[0], e_range[1])

    if function_switch['SPIN_TEXTURE']:
        spin_texture_parameters = INPUT['SPIN_TEXTURE']
        if nspin == 4:
            cal_ST = Spin_Texture(m_tb)
            cal_ST.calculate_spin_texture(**spin_texture_parameters)
        if RANK == 0:
            fermi_energy = input_parameters['fermi_energy']
            cal_ST.print_plot_script(fermi_energy=fermi_energy)

    if function_switch['SURFACE_STATE']:
        surface_state_parameters = INPUT['SURFACE_STATE']
        cal_SS = Surface_State(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_SS.calculate_surface_state(fermi_energy=fermi_energy, **surface_state_parameters)
        if RANK == 0:
            cal_SS.print_plot_script()

    if function_switch['BERRY_CURVATURE']:
        berry_curvature_parameters = INPUT['BERRY_CURVATURE']
        cal_BC = Berry_Curvature(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_BC.calculate_berry_curvature(fermi_energy=fermi_energy, **berry_curvature_parameters)
        if RANK == 0:
            cal_BC.print_plot_script()

    if function_switch['CHERN_NUMBER']:
        chern_number_parameters = INPUT['CHERN_NUMBER']
        cal_CN = Chern_Num(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_CN.calculate_chern_num(fermi_energy=fermi_energy, **chern_number_parameters)

    if function_switch['AHC']:
        ahc_parameters = INPUT['AHC']
        cal_AHC = AHC(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_AHC.calculate_ahc(fermi_energy=fermi_energy, **ahc_parameters)

    if function_switch['ANC']:
        anc_parameters = INPUT['ANC']
        cal_ANC = ANC(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_ANC.calculate_anc(fermi_energy=fermi_energy, **anc_parameters)
        if RANK == 0:
            cal_ANC.print_plot_script()

    if function_switch['OPTICAL_CONDUCTIVITY']:
        optical_conductivity_parameters = INPUT['OPTICAL_CONDUCTIVITY']
        cal_OC = Optical_Conductivity(m_tb)
        cal_OC.calculate_optical_conductivity(**optical_conductivity_parameters)
        if RANK == 0:
            cal_OC.print_plot_script()

    if function_switch['ORBITAL_MAGNETIZATION']:
        orbital_magnetization_parameters = INPUT['ORBITAL_MAGNETIZATION']
        cal_OM = Orbital_Magnetization(m_tb)
        cal_OM.calculate_orbital_magnetization(**orbital_magnetization_parameters)
        if RANK == 0:
            cal_OM.print_plot_script()

    if function_switch['POLARIZATION']:
        polarization_parameters = INPUT['POLARIZATION']
        cal_p = Polarization(m_tb)
        cal_p.calculate_polarization(**polarization_parameters)

    if function_switch['SHIFT_CURRENT']:
        shift_current_parameters = INPUT['SHIFT_CURRENT']
        cal_SC = Shift_Current(m_tb)
        cal_SC.calculate_shift_current(**shift_current_parameters)
        if RANK == 0:
            cal_SC.print_plot_script() 

    if function_switch['WILSON_LOOP']:
        wilson_loop_parameters = INPUT['WILSON_LOOP']
        cal_WL = Wilson_Loop(m_tb)
        cal_WL.calculate_wilson_loop(**wilson_loop_parameters)
        if RANK == 0:
            cal_WL.print_plot_script()

    if function_switch['BERRY_CURVATURE_DIPOLE']:
        berry_curvature_dipole_parameters = INPUT['BERRY_CURVATURE_DIPOLE']
        cal_berry_curvature_dipole = Berry_Curvature_Dipole(m_tb,**input_parameters, **berry_curvature_dipole_parameters)
        cal_berry_curvature_dipole.calculate_berry_curvature_dipole(**berry_curvature_dipole_parameters)
        if RANK == 0:
            cal_berry_curvature_dipole.print_plot_script()
    if function_switch['SHG']:
        shg_parameters = INPUT['SHG']
        cal_shg = Second_Harmonic_Generation(m_tb,**input_parameters, **shg_parameters)
        cal_shg.calculate_shg(**shg_parameters)
        if RANK == 0:
            cal_shg.print_plot_script()
            
    if function_switch['POCKELS']:
        pockels_parameters = INPUT['POCKELS']
        cal_pockels = Pockels(m_tb,**input_parameters, **pockels_parameters)
        cal_pockels.calculate_pockels(**pockels_parameters)
        if RANK == 0:
            cal_pockels.print_plot_script()
        
    if function_switch['CHIRALITY']:
        chirality_parameters = INPUT['CHIRALITY']
        cal_chirality = Chirality(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        cal_chirality.calculate_chirality(fermi_energy=fermi_energy, **chirality_parameters)
        
    if function_switch['CPGE']:
        cpge_parameters = INPUT['CPGE']
        cal_cpge = CPGE(m_tb,**input_parameters, **cpge_parameters)
        cal_cpge.calculate_cpge(**cpge_parameters)
        
    if function_switch['DRUDE_WEIGHT']:
        drude_weight_parameters = INPUT['DRUDE_WEIGHT']
        cal_drude_weight = Drude_Weight(m_tb,**input_parameters, **drude_weight_parameters)
        cal_drude_weight.calculate_drude_weight(**drude_weight_parameters)

    if function_switch['REDUCE_BASIS']:
        reduce_basis_parameters = INPUT['REDUCE_BASIS']
        RBC = Reduce_Basis_Check(m_tb)
        RBC.get_reduce_basis(**reduce_basis_parameters)

    if function_switch['BOLTZ_TRANSPORT']:
        transport_parameters = INPUT['BOLTZ_TRANSPORT']
        TRANS = Boltz_Transport(m_tb)
        fermi_energy = input_parameters['fermi_energy']
        TRANS.calculate_boltz_transport(fermi_energy=fermi_energy, **transport_parameters)
    
    if RANK == 0:
        timer.print_all()
        timer.program_end()


if __name__ == '__main__':
    main()
