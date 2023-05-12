from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.parallel import op_sum
from pyatb.kpt import kpoint_generator
import numpy as np
import time

class grid_integrate_3D:
    def __init__(
        self,
        func,
        origin,
        vect1,
        vect2,
        vect3,
        grid1,
        grid2,
        bar,
        max_point_num = 100
    ):
        self.func = func
        
        self.origin = origin
        self.vect1 = vect1
        self.vect2 = vect2
        self.vect3 = vect3
        self.grid1 = grid1
        self.grid2 = grid2
        self.bar = bar
        self.max_point_num = max_point_num
        
        num_func = func(np.array([origin])).shape
        if num_func == (1,):
            self.num_func = 1
        else:
            self.num_func = num_func[1]

    def integrate(self):
        # first grid sum
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\nEnter the integral solution module ==> \n')
                f.write('\n!!The number of rough k-grid : %d\n'%(self.grid1[0]*self.grid1[1]*self.grid1[2]))

        weight = 1.0 / self.grid1[0] / self.grid1[1] / self.grid1[2]
        ans, point_list = self.__area_sum(self.origin, self.vect1, self.vect2, self.vect3, self.grid1, judge=True)
        ans = ans * weight

        # second adaptive grid sum
        COMM.Barrier()
        if RANK == 0:
            with open(RUNNING_LOG, 'a') as f:
                f.write('\n!!The number of k points for adaptive refinement : %d\n'%(len(point_list)))

        weight = weight / self.grid2[0] / self.grid2[1] / self.grid2[2]
        new_vect1 = self.vect1 / self.grid1[0]
        new_vect2 = self.vect2 / self.grid1[1]
        new_vect3 = self.vect3 / self.grid1[2]
    
        for point in point_list:
            temp_origin = point - (new_vect1 + new_vect2 + new_vect3) * 0.5
            temp_ans = self.__area_sum(temp_origin, new_vect1, new_vect2, new_vect3, self.grid2, judge=False)[0]
            if RANK == 0:
                ans = ans + temp_ans * weight
        
        # result
        self.ans = ans

        COMM.Barrier()
        return ans      

    def __area_sum(self, origin, vect1, vect2, vect3, grid, judge=True):
        point_list = list()
        combine_point_list = list()
        k = kpoint_generator.mp_generator(self.max_point_num, origin, vect1, vect2, vect3, grid)
        ans = 0
        for kpoint in k:
            COMM.Barrier()
            time_start = time.time()

            ik_process = kpoint_generator.kpoints_in_different_process(SIZE, RANK, kpoint)
            kpoint_num = ik_process.k_direct_coor_local.shape[0]
            if kpoint_num:
                ans_list = self.func(ik_process.k_direct_coor_local)
                for i, ele in enumerate(ans_list):
                    if judge:
                        if self.judge_value(ele):
                            point_list.append(ik_process.k_direct_coor_local[i])
                        else:
                            ans = ans + ele
                    else:
                        ans = ans + ele
            
            COMM.Barrier()
            time_end = time.time()
            if RANK == 0:
                with open(RUNNING_LOG, 'a') as f:
                    f.write(' >> Calculated %10d k points, took %.6e s\n'%(kpoint.shape[0], time_end-time_start))

        all_sum = COMM.allreduce(ans, op=op_sum)
        combine_point_list = COMM.allreduce(point_list, op=op_sum)

        return all_sum, combine_point_list

    def judge_value(self, value):
        bar = self.bar

        # Determine if this point needs to be dense, if so, returns true, otherwise return false.
        if self.num_func == 1:
            if abs(value) < bar:
                return False
            else:
                return True
        else:
            if np.linalg.norm(value) < bar:
                return False
            else:
                return True