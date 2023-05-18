import numpy as np
import math
import itertools

class generator:
    def __init__(self,ndim):
        self.ndim = ndim
        
    #generator会根据不同的规则编号与被积函数的维度来产生一组点，以及这组点对应的权重
    #这组点相对中心点的位移有多种
    #编号1 （0，0，0   ，0）型 即中心点 在全空间分布1个
    #编号2 （alpha，，    0）型 在全空间分布有2n个
    #编号3 （beta，beta，，0）型 在全空间分布有2n(n-1)个
    #编号4 （epsilon，epsilon，epsilon，，0）型 在全空间分布有4n(n-1)(n-2)/3个
    
    def rule(self,rule_num):
        self.listlenth = 0
        if rule_num == 5:
            self.rule5()
        elif rule_num == 7:
            self.rule7()
        elif rule_num == 9:
            self.rule9()
        ndim = self.ndim
        totallist = list()
        for j in range(self.wtleng):
            g = self.g[:,j]

            templist = list(set(list(itertools.permutations(g))))
            for i in range(1,ndim):
                temp = [-1]*i
                temp1 = [1]*(ndim-i)
                temp.extend(temp1)
                
                temptemp = list(set(list(itertools.permutations(temp))))
                for templ in temptemp:
                    templ = np.array(templ)
                for templ in temptemp:
                    for k in range(ndim):
                    
                        g[k] = g[k]*templ[k]
                    templist.extend(list(set(list(itertools.permutations(g)))))
            templist = list(set(templist))
            self.rulpts[j] = len(templist)
            
            #print(j,len(templist))
            totallist.extend(templist)
        #print(totallist)
        self.listlenth = len(totallist)
        for i in range(self.listlenth):
            tp = np.array(totallist[i],dtype = float)
            totallist[i] = tp
        self.pointlist = totallist
        


        
    def rule5(self):
        #计算degree5的generator以及weight
        
        ndim = self.ndim 
        wtleng = 6
        self.wtleng = 6
        w = np.zeros([5,wtleng],dtype = float)
        g = np.zeros([ndim,wtleng],dtype = float)
        rulpts = np.zeros(6,dtype = float)
        for i in range(wtleng):
            rulpts[i] = 2*ndim
        
        rulpts[wtleng-1] = 2*ndim*(ndim-1)
        rulpts[0] = 1
        
        lam1 = 3/5
        lam2 = 1/3
        lam3 = 3/4
        lam4 = 4/5
        
        w[0,wtleng-1] = 1/pow((6*lam1),2)
        w[0,wtleng-2] = 1/(6*lam1)-2*(ndim-1)*w[0,wtleng-1]
        
        w[1,wtleng-1] = 1/pow((6*lam1),2)
        w[1,1] = 1/(6*lam3)-2*(ndim-1)*w[1,wtleng-1]*lam1/lam3
        w[2,1] = (1/5-lam4/3)/(2*lam2*(lam2-lam4))
        w[2,3] = (1/5-lam2/3)/(2*lam4*(lam4-lam2))
        w[3,wtleng-1] = 1/(2*ndim*(ndim-1))
        w[4,1] = 1/(4*ndim)
        w[4,3] = 1/(4*ndim)
        
        lam1 = math.sqrt(lam1)
        lam2 = math.sqrt(lam2)
        lam3 = math.sqrt(lam3)
        lam4 = math.sqrt(lam4)
        
        g[0,wtleng-1] = lam1
        g[1,wtleng-1] = lam1
        g[0,wtleng-2] = lam1
        g[0,1] = lam2
        g[0,2] = lam3
        g[0,3] = lam4
        
        for j in range(5):
            w[j,0] = pow(2,ndim)
            for i in range(1,wtleng):
                w[j,i] = pow(2,ndim)*w[j,i]
                w[j,0] = w[j,0] - rulpts[i]*w[j,i]
        self.rulnrm(wtleng,5,rulpts,w)
        self.g = g
        self.w = w
        self.rulpts = rulpts
        return
    
    
    def rule7(self):
        ndim = self.ndim
        wtleng = 6
        self.wtleng = 6
        w = np.zeros([5,wtleng],dtype = float)
        g = np.zeros([ndim,wtleng],dtype = float)
        rulpts = np.zeros(6,dtype = float)
        for i in range(wtleng):
            rulpts[i] = 2*ndim
        twondm = 2**ndim
        rulpts[wtleng-1] = twondm
        rulpts[wtleng-2] = 2*ndim*(ndim-1)
        rulpts[0] = 1
        
        lam0 = 0.4707
        lamp = 0.5625
        lam1 = 4/(15-5/lam0)
        ratio = (1-lam1/lam0)/27
        lam2 = (5-7*lam1-35*ratio)/(7-35*lam1/3-35*ratio/lam0)
        
        w[0,5] = math.pow(1/(3*lam0),3.0)/twondm
        
        
        w[0,5] = 1/ (3*lam0)**3/twondm
        w[0,4] = (1-5*lam0/3)/ (60* (lam1-lam0)*lam1**2)
        w[0,2] = (1-5*lam2/3-5*twondm*w[0,5]*lam0* (lam0-lam2))/(10*lam1* (lam1-lam2)) - 2* (ndim-1)*w[0,4]
        w[0,1] = (1-5*lam1/3-5*twondm*w[0,5]*lam0* (lam0-lam1))/(10*lam2* (lam2-lam1))

        w[1,5] = 1/ (36*lam0**3)/twondm
        w[1,4] = (1-9*twondm*w[1,5]*lam0**2)/ (36*lam1**2)
        w[1,2] = (1-5*lam2/3-5*twondm*w[1,5]*lam0* (lam0-lam2))/(10*lam1* (lam1-lam2)) - 2* (ndim-1)*w[1,4]
        w[1,1] = (1-5*lam1/3-5*twondm*w[1,5]*lam0* (lam0-lam1))/(10*lam2* (lam2-lam1))
        w[2,5] = 5/ (108*lam0**3)/twondm
        w[2,4] = (1-9*twondm*w[2,5]*lam0**2)/ (36*lam1**2)
        w[2,2] = (1-5*lamp/3-5*twondm*w[2,5]*lam0* (lam0-lamp))/(10*lam1* (lam1-lamp)) - 2* (ndim-1)*w[2,4]
        w[2,3] = (1-5*lam1/3-5*twondm*w[2,5]*lam0* (lam0-lam1))/(10*lamp* (lamp-lam1))
        w[3,5] = 1/ (54*lam0**3)/twondm
        w[3,4] = (1-18*twondm*w[3,5]*lam0**2)/ (72*lam1**2)
        w[3,2] = (1-10*lam2/3-10*twondm*w[3,5]*lam0* (lam0-lam2))/(20*lam1* (lam1-lam2)) - 2* (ndim-1)*w[3,4]
        w[3,1] = (1-10*lam1/3-10*twondm*w[3,5]*lam0* (lam0-lam1))/(20*lam2* (lam2-lam1))

        lam0 = math.sqrt(lam0)
        lam1 = math.sqrt(lam1)
        lam2 = math.sqrt(lam2)
        lamp = math.sqrt(lamp)
        
        
        for i in range(ndim):
            g[i,wtleng-1] = lam0
          
     
        g[0,wtleng-2] = lam1
        g[1,wtleng-2] = lam1
        g[0,wtleng-5] = lam2
        g[0,wtleng-4] = lam1
        g[0,wtleng-3] = lamp
        for j in range(5):
            w[j,0] = twondm
            for i in range(1,wtleng):
                w[j,i] = twondm*w[j,i]
                w[j,0] = w[j,0]-rulpts[i]*w[j,i]
        self.rulnrm ( wtleng, 5, rulpts, w )
        
        self.g = g
        self.w = w
        self.rulpts = rulpts
        return
    
    
    def rule9(self):
        wtleng = 9
        self.wtleng = 9
        ndim = self.ndim
        if ndim==2 :
            wtleng = 8
        w = np.zeros([5,wtleng],dtype = float)
        g = np.zeros([ndim,wtleng],dtype = float)
        rulpts = np.zeros(wtleng,dtype = float)
        for j in range(wtleng):
            rulpts[j] = 2*ndim
        twondm = 2**ndim
        rulpts[wtleng-1] = twondm
        if ndim>2:
            rulpts[7] = (4*ndim* (ndim-1)* (ndim-2))/3
        rulpts[7] = 4*ndim* (ndim-1)
        rulpts[5] = 2*ndim* (ndim-1)
        rulpts[0] = 1

        lam0 = 0.4707
        lam1 = 4/ (15-5/lam0)
        ratio = (1-lam1/lam0)/27
        lam2 = (5-7*lam1-35*ratio)/ (7-35*lam1/3-35*ratio/lam0)
        ratio = ratio* (1-lam2/lam0)/3
        lam3 = (7-9* (lam2+lam1)+63*lam2*lam1/5-63*ratio)/(9-63* (lam2+lam1)/5+21*lam2*lam1-63*ratio/lam0)
        lamp = 0.0625

        
        w[0,wtleng-1] = 1/ math.pow((3*lam0),4)/twondm
        
        if ndim>2:
            w[0,7] = (1-1/ (3*lam0))/ math.pow((6*lam1),3)
     
        w[0,6] = (1-7* (lam0+lam1)/5+7*lam0*lam1/3)/(84*lam1*lam2* (lam2-lam0)* (lam2-lam1))
        w[0,5] = (1-7* (lam0+lam2)/5+7*lam0*lam2/3)/(84*lam1*lam1* (lam1-lam0)* (lam1-lam2)) -w[0,6]*lam2/lam1 - 2* (ndim-2)*w[0,7]
        w[0,3] = (1-9* ((lam0+lam1+lam2)/7- (lam0*lam1+lam0*lam2+lam1*lam2)/5)-3*lam0*lam1*lam2)/(18*lam3* (lam3-lam0)* (lam3-lam1)* (lam3-lam2))
        w[0,2] = (1-9* ((lam0+lam1+lam3)/7- (lam0*lam1+lam0*lam3+lam1*lam3)/5)-3*lam0*lam1*lam3)/(18*lam2* (lam2-lam0)* (lam2-lam1)* (lam2-lam3)) -2* (ndim-1)*w[0,6]
        w[0,1] = (1-9* ((lam0+lam2+lam3)/7- (lam0*lam2+lam0*lam3+lam2*lam3)/5)-3*lam0*lam2*lam3)/(18*lam1* (lam1-lam0)* (lam1-lam2)* (lam1-lam3)) -2* (ndim-1)* (w[0,6]+w[0,5]+ (ndim-2)*w[0,7])

        w[1,wtleng-1] = 1/ (108*lam0**4)/twondm
        
        if(ndim>2):
            w[1,7]= (1-27*twondm*w[1,8]*lam0**3)/ (6*lam1)**3
      
        w[1,6] = (1-5*lam1/3-15*twondm*w[1,wtleng-1]*lam0**2* (lam0-lam1))/(60*lam1*lam2* (lam2-lam1))
        w[1,5] = (1-9* (8*lam1*lam2*w[1,6]+twondm*w[1,wtleng-1]*lam0**2))/(36*lam1*lam1) - 2*w[1,7]* (ndim-2)
        w[1,3] = (1-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w[1,wtleng-1]*lam0* (lam0-lam1)* (lam0-lam2)))/(14*lam3* (lam3-lam1)* (lam3-lam2))
        w[1,2] = (1-7* ((lam1+lam3)/5-lam1*lam3/3+twondm*w[1,wtleng-1]*lam0* (lam0-lam1)* (lam0-lam3)))/(14*lam2* (lam2-lam1)* (lam2-lam3)) - 2* (ndim-1)*w[1,6]
        w[1,1] = (1-7* ((lam2+lam3)/5-lam2*lam3/3+twondm*w[1,wtleng-1]*lam0* (lam0-lam2)* (lam0-lam3)))/(14*lam1* (lam1-lam2)* (lam1-lam3)) -2* (ndim-1)* (w[1,6]+w[1,5]+ (ndim-2)*w[1,7])
        w[2,wtleng-1] = 5/ (324*lam0**4)/twondm
        
        if ndim>2:
            w[2,7] = (1-27*twondm*w[2,8]*lam0**3)/ (6*lam1)**3
      
        w[2,6] = (1-5*lam1/3-15*twondm*w[2,wtleng-1]*lam0**2* (lam0-lam1))/(60*lam1*lam2* (lam2-lam1))
        w[2,5] = (1-9* (8*lam1*lam2*w[2,6]+twondm*w[2,wtleng-1]*lam0**2))/(36*lam1*lam1) - 2*w[2,7]* (ndim-2)
        w[2,4] = (1-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w[2,wtleng-1]*lam0* (lam0-lam1)* (lam0-lam2)))/(14*lamp* (lamp-lam1)* (lamp-lam2))
        w[2,2] = (1-7* ((lam1+lamp)/5-lam1*lamp/3+twondm*w[2,wtleng-1]*lam0* (lam0-lam1)* (lam0-lamp)))/(14*lam2* (lam2-lam1)* (lam2-lamp)) - 2* (ndim-1)*w[2,6]
        w[2,1] = (1-7* ((lam2+lamp)/5-lam2*lamp/3+twondm*w[2,wtleng-1]*lam0* (lam0-lam2)* (lam0-lamp)))/(14*lam1* (lam1-lam2)* (lam1-lamp)) -2* (ndim-1)* (w[2,6]+w[2,5]+ (ndim-2)*w[2,7])
        w[3,wtleng-1] = 2/ (81*lam0**4)/twondm
        
        if ndim>2:
            w[3,7] = (2-27*twondm*w[3,8]*lam0**3)/ (6*lam1)**3
      
        w[3,6] = (2-15*lam1/9-15*twondm*w[3,wtleng-1]*lam0* (lam0-lam1))/(60*lam1*lam2* (lam2-lam1))
        w[3,5] = (1-9* (8*lam1*lam2*w[3,6]+twondm*w[3,wtleng-1]*lam0**2))/(36*lam1*lam1) - 2*w[3,7]* (ndim-2)
        w[3,3] = (2-7* ((lam1+lam2)/5-lam1*lam2/3+twondm*w[3,wtleng-1]*lam0* (lam0-lam1)* (lam0-lam2)))/(14*lam3* (lam3-lam1)* (lam3-lam2))
        w[3,2] = (2-7* ((lam1+lam3)/5-lam1*lam3/3+twondm*w[3,wtleng-1]*lam0* (lam0-lam1)* (lam0-lam3)))/(14*lam2* (lam2-lam1)* (lam2-lam3)) - 2* (ndim-1)*w[3,6]
        w[3,1] = (2-7* ((lam2+lam3)/5-lam2*lam3/3+twondm*w[3,wtleng-1]*lam0* (lam0-lam2)* (lam0-lam3)))/(14*lam1* (lam1-lam2)* (lam1-lam3)) -2* (ndim-1)* (w[3,6]+w[3,5]+ (ndim-2)*w[3,7])
        w[4,1] = 1/ (6*lam1)
        
        lam0 = math.sqrt(lam0)
        lam1 = math.sqrt(lam1)
        lam2 = math.sqrt(lam2)
        lam3 = math.sqrt(lam3)
        lamp = math.sqrt(lamp)
        
        for i in range(ndim):
            g[i,wtleng-1] = lam0
        if ndim>2:
            g[0,7] = lam1
            g[1,7] = lam1
            g[2,7] = lam1
      
        g[0,6] = lam1
        g[1,6] = lam2
        g[0,5] = lam1
        g[1,5] = lam1
        g[0,4] = lamp
        g[0,3] = lam3
        g[0,2] = lam2
        g[0,1] = lam1
        
        for j in range(5):
            w[j,0] = twondm
            for i in range(1,wtleng):
                w[j,i] = twondm*w[j,i]
                w[j,0] = w[j,0] - rulpts[i]*w[j,i]
      
        self.rulnrm ( wtleng, 5, rulpts, w )
        self.g = g
        self.w = w
        self.rulpts = rulpts
        return

    
    
    
    
    
    def rulnrm(self,lenrul,numnul,rulpts,w):
        normcf = 0
        normnl = 0
        for i in range(lenrul):
            normcf = normcf+rulpts[i]*w[0,i]*w[0,i]
        for k in range(1,numnul):
            for i in range(lenrul):
                w[k,i] = w[k,i]-w[0,i]
            for j in range(1,k-1):
                alpha = 0
                for i in range(lenrul):
                    alpha = alpha+rulpts[i]*w[j,i]*w[k,i]
                for i in range(lenrul):
                    w[k,i] = w[k,i]+alpha*w[j,i]
            for i in range(lenrul):
                normnl = normnl+rulpts[i]*w[k,i]*w[k,i]
            alpha = math.sqrt(normcf/normnl)
            for i in range(lenrul):
                w[k,i] = alpha*w[k,i]
        self.w = w
        return


