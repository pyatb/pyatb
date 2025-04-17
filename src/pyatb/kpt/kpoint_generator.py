import numpy as np

class kpoints_in_different_process:
    def __init__(self, process_num, rank_id, k_direct_coor_global):
        self.process_num = process_num
        self.rank_id = rank_id
        self.k_direct_coor_global = k_direct_coor_global
        self.k_direct_coor_local = None
        self.ik_start_index = 0
        self.__distribute_kpoint_to_different_process()

    def __distribute_kpoint_to_different_process(self):
        kpoint_num_global = self.k_direct_coor_global.shape[0]
        remain = kpoint_num_global % self.process_num
        kpoint_num = kpoint_num_global // self.process_num
        if self.rank_id < remain:
            kpoint_num = kpoint_num + 1
        self.k_direct_coor_local = np.empty([kpoint_num, 3], dtype=float)

        if self.rank_id < remain:
            self.ik_start_index = self.rank_id * kpoint_num
        else:
            self.ik_start_index = self.rank_id * kpoint_num + remain

        jk = self.ik_start_index
        for ik in range(kpoint_num):
            self.k_direct_coor_local[ik] = self.k_direct_coor_global[jk]
            jk = jk + 1


class string_in_different_process:
    def __init__(self, process_num, rank_id, string_direct_coor_global):
        self.process_num = process_num
        self.rank_id = rank_id
        self.string_direct_coor_global = string_direct_coor_global
        self.string_direct_coor_local = None
        self.start_index = 0

        self.__distribute_kpoint_to_different_process()

    def __distribute_kpoint_to_different_process(self):
        tem_shape2 = self.string_direct_coor_global.shape[1]
        string_num_global = self.string_direct_coor_global.shape[0]
        remain = string_num_global % self.process_num
        string_num_local = string_num_global // self.process_num
        if self.rank_id < remain:
            string_num_local = string_num_local + 1
        self.string_direct_coor_local = np.empty([string_num_local, tem_shape2, 3], dtype=float)

        if self.rank_id < remain:
            self.start_index = self.rank_id * string_num_local
        else:
            self.start_index = self.rank_id * string_num_local + remain

        jk = self.start_index
        for ik in range(string_num_local):
            self.string_direct_coor_local[ik] = self.string_direct_coor_global[jk]
            jk = jk + 1


class mp_generator:
    def __init__(self, max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid):
        self.max_kpoint_num = max_kpoint_num
        self.k_start = k_start
        self.k_vect1 = k_vect1
        self.k_vect2 = k_vect2
        self.k_vect3 = k_vect3
        self.grid = grid
        self.total_kpoint_num = grid[0] * grid[1] * grid[2]
        self.current_kpoint_num = 0

    def __get_kpoint(self):
        grid = self.grid
        remain_kpoint_num = self.total_kpoint_num - self.current_kpoint_num
        if remain_kpoint_num > self.max_kpoint_num:
            k_direct_coor = np.empty([self.max_kpoint_num, 3], dtype=float)
        else:
            k_direct_coor = np.empty([remain_kpoint_num, 3], dtype=float)

        total_count = 0
        count = 0
        for x in range(grid[0]):
            for y in range(grid[1]):
                for z in range(grid[2]):
                    total_count = total_count + 1
                    if total_count > self.current_kpoint_num:
                        if count == self.max_kpoint_num:
                            self.current_kpoint_num = self.current_kpoint_num + count
                            return k_direct_coor
                        k_direct_coor[count] = self.k_start + self.k_vect1 * (x/grid[0]) + self.k_vect2 * (y/grid[1]) + self.k_vect3 * (z/grid[2])
                        count = count + 1
        self.current_kpoint_num = self.current_kpoint_num + count
        return k_direct_coor

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_kpoint_num < self.total_kpoint_num: 
            return self.__get_kpoint()
        else:
            raise StopIteration

class line_generator:
    def __init__(self, max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line):
        self.max_kpoint_num = max_kpoint_num
        self.high_symmetry_kpoint = high_symmetry_kpoint
        self.kpoint_num_in_line = kpoint_num_in_line
        self.total_kpoint_num = kpoint_num_in_line.sum() - kpoint_num_in_line[-1] + 1
        self.current_kpoint_num = 0

    def __get_kpoint(self):
        remain_kpoint_num = self.total_kpoint_num - self.current_kpoint_num
        if remain_kpoint_num > self.max_kpoint_num:
            k_direct_coor = np.empty([self.max_kpoint_num, 3], dtype=float)
        else:
            k_direct_coor = np.empty([remain_kpoint_num, 3], dtype=float)

        total_count = 0
        count = 0
        for i in range(1, self.high_symmetry_kpoint.shape[0]):
            dk = (self.high_symmetry_kpoint[i] - self.high_symmetry_kpoint[i-1]) / self.kpoint_num_in_line[i-1]
            for j in range(self.kpoint_num_in_line[i-1]):
                total_count = total_count + 1
                if total_count > self.current_kpoint_num:
                    k_direct_coor[count] = self.high_symmetry_kpoint[i-1] + dk * j
                    count = count + 1
                    if count == self.max_kpoint_num:
                        self.current_kpoint_num = self.current_kpoint_num + count
                        return k_direct_coor
        k_direct_coor[count] = self.high_symmetry_kpoint[-1]
        self.current_kpoint_num = self.current_kpoint_num + count + 1
        return k_direct_coor

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_kpoint_num < self.total_kpoint_num: 
            return self.__get_kpoint()
        else:
            raise StopIteration

class array_generater:
    def __init__(self, max_kpoint_num, total_kpoint):
        self.max_kpoint_num = max_kpoint_num
        
        self.total_kpoint = total_kpoint
        self.total_kpoint_num = self.total_kpoint.shape[0]
        self.current_kpoint_num = 0

    def __get_kpoint(self):
        remain_kpoint_num = self.total_kpoint_num - self.current_kpoint_num
        if remain_kpoint_num > self.max_kpoint_num:
            kpoint_num = self.max_kpoint_num
            k_direct_coor = np.empty([self.max_kpoint_num, 3], dtype=float)
        else:
            kpoint_num = remain_kpoint_num
            k_direct_coor = np.empty([remain_kpoint_num, 3], dtype=float)
        
        k_direct_coor = self.total_kpoint[self.current_kpoint_num:self.current_kpoint_num+kpoint_num, :]
         
        self.current_kpoint_num = self.current_kpoint_num + kpoint_num
        return k_direct_coor

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_kpoint_num < self.total_kpoint_num: 
            return self.__get_kpoint()
        else:
            raise StopIteration


class string_generator:
    def __init__(self, max_kpoint_num, k_start, k_vect1, k_vect2, nk1, nk2):
        self.max_kpoint_num = max_kpoint_num
        self.k_start = k_start
        self.k_vect1 = k_vect1
        self.k_vect2 = k_vect2
        self.nk1 = nk1  # the number of k points on each string
        self.nk2 = nk2  # total string number
        self.string_num = max_kpoint_num // nk1
        self.current_string_num = 0

    def __get_string(self):
        remain = self.nk2 - self.current_string_num
        if remain > self.string_num:
            tem_string_num = self.string_num
            string_direct_coor = np.zeros([self.string_num, self.nk1, 3], dtype=float)
        else:
            tem_string_num = remain
            string_direct_coor = np.zeros([remain, self.nk1, 3], dtype=float)
        
        for ik2 in range(self.current_string_num, self.current_string_num+tem_string_num):
            for ik1 in range(self.nk1):
                string_direct_coor[ik2-self.current_string_num, ik1] = self.k_start + self.k_vect1 * ik1 / self.nk1 + self.k_vect2 * ik2 / (self.nk2 - 1)

        self.current_string_num = self.current_string_num + tem_string_num

        return string_direct_coor

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_string_num < self.nk2:
            return self.__get_string()
        else:
            raise StopIteration

class string_generator_3d:
    def __init__(self, max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, nk1, nk2, nk3, direction):
        self.max_kpoint_num = max_kpoint_num
        self.k_start = k_start

        if direction == 0:
            self.k_vect1 = k_vect2
            self.k_vect2 = k_vect3
            self.k_vect3 = k_vect1
            self.grid = [nk2, nk3, nk1]
        elif direction == 1:
            self.k_vect1 = k_vect1
            self.k_vect2 = k_vect3
            self.k_vect3 = k_vect2
            self.grid = [nk1, nk3, nk2]
        elif direction == 2:
            self.k_vect1 = k_vect1
            self.k_vect2 = k_vect2
            self.k_vect3 = k_vect3
            self.grid = [nk1, nk2, nk3]

        self.direction = direction
        self.total_string_num = self.grid[0] * self.grid[1]
        self.string_num = max_kpoint_num // self.grid[2]
        self.current_string_num = 0

    def __get_string(self):
        remain = self.total_string_num - self.current_string_num
        if remain > self.string_num:
            tem_string_num = self.string_num
            string_direct_coor = np.zeros([self.string_num, self.grid[2], 3], dtype=float)
        else:
            tem_string_num = remain
            string_direct_coor = np.zeros([remain, self.grid[2], 3], dtype=float)
        
        count = 0
        for ik1 in range(self.grid[0]):
            for ik2 in range(self.grid[1]):
                if count >= self.current_string_num and count < self.current_string_num + tem_string_num:
                    for ik3 in range(self.grid[2]):
                        string_direct_coor[count-self.current_string_num, ik3] = (
                            self.k_start + self.k_vect1 * ik1 / self.grid[0] + self.k_vect2 * ik2 / self.grid[1] + self.k_vect3 * ik3 / self.grid[2]
                        )
                count = count + 1
        

        self.current_string_num = self.current_string_num + tem_string_num

        return string_direct_coor

    def __iter__(self):
        return self

    def __next__(self):
        if self.current_string_num < self.total_string_num:
            return self.__get_string()
        else:
            raise StopIteration



'''
k_start = np.array([0, 0, 0], dtype=float)
k_vect1 = np.array([1, 0, 0], dtype=float)
k_vect2 = np.array([0, 1, 0], dtype=float)
k_vect3 = np.array([0, 0, 1], dtype=float)
grid = np.array([4, 4, 4], dtype=int)
max_kpoint_num = 10

high_symmetry_kpoint = np.array([[0, 0, 0], [1, 1, 1]], dtype=float)
kpoint_num_in_line = np.array([20, 1], dtype=int)


k = mp_generator(max_kpoint_num, k_start, k_vect1, k_vect2, k_vect3, grid)

# k = line_generator(max_kpoint_num, high_symmetry_kpoint, kpoint_num_in_line)

for kpoint in k:
    ik_process_0 = kpoints_in_different_process(2, 0, kpoint)
    ik_process_1 = kpoints_in_different_process(2, 1, kpoint)
    print(ik_process_0.k_direct_coor_local)
    print(ik_process_1.k_direct_coor_local)
'''