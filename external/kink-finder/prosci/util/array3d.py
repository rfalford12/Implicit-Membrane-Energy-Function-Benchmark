#
#  array3d.py
#  
#
#  Created by Sebastian Kelm on 18/09/2010.
#  Copyright (c) 2010 Sebastian Kelm. All rights reserved.
#

class Array3D:
    def __init__(self, shape):
        "Linear implementation of a 3d array (of lists)"
        assert len(shape) == 3
        self.shape = shape
        self.s0 = shape[1]*shape[2]
        self.s1 = shape[1]
        assert isinstance(self.s0, int)
        assert isinstance(self.s1, int)
        self.data = [None] * (self.adress2index(self.shape-1)+1)
    
    def adress2index(self, a):
        return a[0]*self.s0 + a[1]*self.s1 + a[2]
    
    def __getitem__(self, a):
        return self.data[self.adress2index(a)]

    #def __setitem__(self, a, val):
    #    return self.data[self.adress2index(a)] = val
    
    def add(self, a, val):
        """Array3D.add(address, value) : adds the value to the given address in the array.
        
        For each array address, values are saved in a list, to allow multiple values.
        """
        i = self.adress2index(a)
        L = self.data[i]
        if L is None:
            self.data[i] = [val]
        else:
            L.append(val)
