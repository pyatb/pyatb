import operator
import numpy as np
from scipy.sparse import csc_matrix


class multiXR:
    """
    Class for holding XR (eg HR, SR, rR) matrices.

    Attributes
    ----------
    XR_des : tuple
        Describe the type of multiXR.

    des : str
        One of 'H', 'S', 'r_x', 'r_y', 'r_z'.

    R_num : int
        The number of R indicators.

    R_direct_coor : np.ndarray[int]
        shape=(R_num, 3). R's fractional coordinates.

    basis_num : int
        Matrix dimension of XR.

    XR : scipy.sparse.csc_matrix
        shape=(R_num, (basis_num + 1)*basis_num / 2). Data of HR or SR or rR. Each row of the sparse matrix corresponds to the upper triangular 
        part of each X[R]. XR's unit is eV for HR, is angstrom for rR.
    """

    XR_des = ('H', 'S', 'r_x', 'r_y', 'r_z')

    def __init__(self, des):
        # Is it HR, SR, or rR
        if des in multiXR.XR_des:
            self.des = des
        else:
            raise ValueError("multiXR des must be one of 'H', 'S', 'r_x', 'r_y', 'r_z'")
        
    def set_XR(self, R_num, R_direct_coor, basis_num, XR):
        """
        Set the data related to XR, the number and fractional coordinates of R, 
        the dimension of the matrix (square matrix) corresponding to each R, 
        and the specific value of XR.

        Parameters
        ----------
        R_num : int, 
            the number of R indicators.

        R_direct_coor : np.ndarray[int]
            shape=(R_num, 3). R's fractional coordinates.

        basis_num : int 
            matrix dimension of XR.

        XR : scipy.sparse.csc_matrix[complex] or np.ndarray[complex]
            if the type is csc_matrix, shape=(R_num, (basis_num + 1)*basis_num / 2); if the type is 
            ndarray, shape=(R_num, basis_num, basis_num). Data of HR or SR or rR. Each row of the sparse 
            matrix corresponds to the upper triangular part of each X[R]. XR's unit is eV for HR, 
            is angstrom for rR.
        """
        self.R_num = R_num
        self.R_direct_coor = R_direct_coor
        self.basis_num = basis_num
        
        XR_shape = (R_num, int((basis_num + 1)*basis_num / 2))
        if isinstance(XR, csc_matrix) and operator.eq(XR.shape, XR_shape):
            self.XR = XR
        elif isinstance(XR, np.ndarray) and operator.eq(XR.shape, (R_num, basis_num, basis_num)):
            temp_XR = np.zeros(XR_shape, dtype=complex)
            for iR in range(R_num):
                count = 0
                for row in range(basis_num):
                    for col in range(row, basis_num):
                        temp_XR[iR, count] = XR[row, col]
                        count = count + 1
            self.XR = csc_matrix(temp_XR)
