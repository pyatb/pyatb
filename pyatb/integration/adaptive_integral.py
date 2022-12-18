from pyatb.integration.generator import generator
from pyatb.kpt import kpoint_generator
from pyatb import RANK, COMM, SIZE, OUTPUT_PATH, RUNNING_LOG, timer
from pyatb.parallel import op_gather_numpy


import numpy as np
import math
import time
import os

class adaptive_integral:
    def __init__(
        self,
        func,
        start,
        end,
        initial_slice=np.array([1,1,1],dtype = int),
        eps_abs = 1e-3,
        eps_rel = 1e-2,
        numfun = 1,
        max_kpoint_num = 5000,
        buffer_max = 400,
        output_path = os.getcwd() + '/'+'integraton.log'):
        
        self.gen = generator(3)
        self.gen.rule(7)
        
        self.start = start
        self.end = end
        
        self.output_path = output_path
        
        self.func = func
        self.numfun = numfun
        
        self.subregion_list = list()
        self.initial_slice = initial_slice
        self.waiting_list = list()
        self.init_slice(self.initial_slice)
        
        self.max_kpoint_num = max_kpoint_num
        
        self.eps_abs = eps_abs
        self.eps_rel = eps_rel
        
        self.buffer_max = buffer_max
        #最多同时计算buffer_max个区域
        self.step_max = 500
        #最多运行的步数
        self.answer_list = list()
        #统计所有的答案
        self.ndim = 3
        return
    def init_slice(self,grid):
        delta = self.end-self.start
        v = np.zeros(3,dtype = float)
        for i in range(3):
            v[i] = delta[i]/grid[i]
        for i in range(grid[0]):
            for j in range(grid[1]):
                for k in range(grid[2]):
                    temp_start = self.start+np.array([i*v[0],j*v[1],k*v[2]])
                    temp_end = temp_start+v
                    
                    self.waiting_list.append(self.subregion(temp_start,temp_end))
        return
    def integrate(self):
        #整体计算的函数
        #首先将整体进行切分 将切分后的区域加入waitinglist
        #运行regionalintegrate计算waitinglist中的点
        #将subregionlist中前一部分的区域进行切分
        if RANK==0:
            with open(self.output_path,'w') as f:
                print("use self adaptive integration method",file = f)
                localtime = time.asctime(time.localtime(time.time()))
                print("start integration     "+str(localtime),file = f)
        for i in range(self.step_max):
            if RANK == 0:
                with open(self.output_path,'a') as f:
                    localtime = time.asctime(time.localtime(time.time()))
                    print(str(i)+' generation    '+str(localtime),file = f)
                    print('integration region:'+str(len(self.waiting_list)),file = f)
            self.regional_integrate()
            temp_err = 0
            [err,val] = self.cal_err_val()#计算总的值和误差
            self.ans = val
            self.answer_list.append(self.ans)
            if len(self.answer_list)>2:
                temp_err = np.fabs((self.answer_list[-2]-self.answer_list[-1])/2)
            total_err = np.maximum(err,temp_err)
            if RANK== 0:
                
                with open(self.output_path,'a') as f:
                    print(total_err,val,file = f)

            if np.linalg.norm(total_err) <= self.eps_abs and np.linalg.norm(total_err)<=self.eps_rel*np.linalg.norm(val):
                self.ans = val
                with open(self.output_path,'a') as f:
                    localtime = time.asctime(time.localtime(time.time()))
                    print('integration complete!    '+str(localtime),file = f)
                    print(total_err,val,file = f)
                return val
            work_lenth = 0
            for region in self.subregion_list:
                temp_err = region.err/region.V*(self.end[0]-self.start[0])*(self.end[1]-self.start[1])*(self.end[2]-self.start[2])
                if np.linalg.norm(temp_err) > self.eps_abs and np.linalg.norm(err):
                    work_lenth=work_lenth+1
            work_lenth = min(self.buffer_max,work_lenth)
            if work_lenth == 0:
                self.ans = val
                with open(self.output_path,'a') as f:
                    localtime = time.asctime(time.localtime(time.time()))
                    print('integration complete!    '+str(localtime),file = f)
                    print(total_err,val,file = f)
                return val
            self.waiting_list.clear()
            temp_region_list = self.subregion_list[:work_lenth].copy()
            del self.subregion_list[:work_lenth]
            temp_waiting_list = list()
            for region in temp_region_list:
                div = region.div.copy()
                direction = 0
                dmax = div[direction]
                for i in range(self.ndim):
                    if dmax<div[i]:
                        dmax = div[i]
                        direction = i
                
                start1 = region.start.copy()
                end2 = region.end.copy()
                width = region.width.copy()
                width = width.astype(float)
                #print(1,width[direction])
                width[direction]=width[direction]*0.5
                #print(2,width[direction])
                end1 = start1+width
                start2 = end2-width
                temp_waiting_list.append(self.subregion(start1,end1))
                temp_waiting_list.append(self.subregion(start2,end2))
            self.waiting_list = temp_waiting_list
    def cal_err_val(self):
        err = 0
        val = 0
        for region in self.subregion_list:
            #print(region.start,region.end)
            err = err+region.err
            val = val+region.basval
        
        return[err,val]
    def regional_integrate(self):
        pointlist = list()
        #临时存储需要计算的点
        for region in self.waiting_list:
            
            for g in self.gen.pointlist:
                point = np.zeros(self.ndim,dtype = float)
                for i in range(self.ndim):
                    point[i] = region.centre[i]+g[i]*region.hwidth[i]
                pointlist.append(point)
                
        #将待计算区域中的所有点都加入pointlist当中
        kpoint_total = np.array(pointlist)
        ans_list_all = None
        k = kpoint_generator.array_generater(self.max_kpoint_num,kpoint_total)
        for kpoint in k:
            ik_process = kpoint_generator.kpoints_in_different_process(SIZE,RANK, kpoint)
            ans_list = self.func(ik_process.k_direct_coor_local)
            tem_ans_list = COMM.reduce(ans_list, root=0, op=op_gather_numpy)
            
            if RANK == 0:
                if ans_list_all is None:
                    ans_list_all = tem_ans_list
                else:
                    ans_list_all = np.r_[ans_list_all, tem_ans_list]

        value_list = COMM.bcast(ans_list_all, root=0)
        
        #传递给func计算这些点的值
        
        startpoint = 0
        endpoint = self.gen.listlenth
        #使用startpoint和endpoint标记该区域中的计算过的第一个点和最后一个点
        #然后对每个区域做循环，计算区域中的各种值
        
        for i in range(len(self.waiting_list)):
            [great,basval,div] = self.evaluate(self.waiting_list[i],value_list[startpoint:endpoint])
            self.waiting_list[i].err = great
            self.waiting_list[i].basval = basval
            self.waiting_list[i].div = div
            if (len(self.subregion_list))==0:
                self.subregion_list.append(self.waiting_list[i])
            else:
                for j in range(len(self.subregion_list)):
                    if self.subregion_list[j].err<great:
                        self.subregion_list.insert(j,self.waiting_list[i])
                        break
                    elif j==(len(self.subregion_list)-1):
                        self.subregion_list.append(self.waiting_list[i])
        
            startpoint = startpoint+self.gen.listlenth
            endpoint = endpoint+self.gen.listlenth  
        return
    def evaluate(self,region,valuelist):
        vol = 1.0
        ndim = self.ndim
        w = self.gen.w
        #print(region.hwidth)
        for i in range(ndim):
            vol = vol*region.hwidth[i]
        wtleng = self.gen.wtleng
        numfun = self.numfun
        rgnerr = np.zeros([numfun,wtleng])
        num = 0
        numnul = 4
        basval = np.zeros([numfun])
        null = np.zeros([numfun,numnul])
        for i in range(wtleng):
            for k in range(int(self.gen.rulpts[i])):
                rgnerr[:,i] = rgnerr[:,i]+valuelist[num]
                num = num+1
        for i in range(wtleng):
            basval = basval+self.gen.w[0,i]*rgnerr[:,i]
            for k in range(numnul):
                null[:,k] = null[:,k]+w[k+1,i]*rgnerr[:,i]
        
            great = 0
        basval= basval*vol
        #print(vol)

        rgnerr_total = np.zeros(numfun)
        for j in range(numfun):
            rgnerr_total[j] = self.twonrm(null[j,0],null[j,1])
            rgncmp = self.twonrm(null[j,1],null[j,2])
            rgncpt = self.twonrm(null[j,2],null[j,3])

            if(4*rgnerr_total[j]<rgncmp and 2*rgncmp < rgncpt):
                rgnerr_total[j] = rgnerr_total[j]/2
            if(2*rgnerr_total[j]>rgncmp):
                rgnerr_total[j] = max(rgnerr_total[j],rgncmp)
            rgnerr_total[j] = vol*rgnerr_total[j]*(i+abs(basval[j]))
            great = great+rgnerr_total[j]
        great = great*vol
      
        r1 = max(np.fabs(self.gen.pointlist[1]))
        r2 = max(np.fabs(self.gen.pointlist[int(1+self.gen.rulpts[1])]))
        v3 = np.sum(np.abs(valuelist[0]))*np.ones([self.ndim])
        v2 = np.zeros([self.ndim])
        v1 = np.zeros([self.ndim])
        for i in range(1,int(1+self.gen.rulpts[1])):
            v1 = v1+np.sum(np.abs(valuelist[i]))*np.abs(self.gen.pointlist[i])
        v1= v1/r1
        for i in range(int(1+self.gen.rulpts[1]),int(1+self.gen.rulpts[2]+self.gen.rulpts[1])):
            v2 = v2+np.sum(np.abs(valuelist[i]))*np.abs(self.gen.pointlist[i])
        v2 = v2/r2
        div = v1/(r1*r1)-v2/(r2*r2)+(1/(r2*r2)-1/(r1*r1))*v3*2
        
        div = np.abs(div)
        #print(great,basval,div)
        return [great,basval,div]
    def twonrm(self,x,y):
        sqttwo = 1.4142135623730950488
        abx = abs(x)
        aby = abs(y)
        if abx>aby:
            twonrm = abx*math.sqrt(1+(aby/abx)**2)
        elif abx<aby:
            twonrm = aby*math.sqrt(1+(abx/aby)**2)
        else:
            twonrm = abx*sqttwo
        return twonrm

    class subregion:
        #存储这个区域的基本信息
        #积分上下限，区域中心
        #积分值、积分误差
        def __init__(self,start,end):
            self.start = start
            self.end = end
            self.centre = (start+end)/2.0
            self.width = end-start
            self.hwidth = self.width/2.0
            self.V = self.width[0]*self.width[1]*self.width[2]
            return
        