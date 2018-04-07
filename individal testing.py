import numpy as np

def absangle(pt1x,pt1y,pt2x,pt2y):
    angle = []
    for i in range(0,len(pt1x)):
        angle.append(np.arctan2((pt2y[i]-pt1y[i]),(pt2x[i]-pt1x[i])))
    return angle

pt1x = [5,6,9,4,5]
pt2x = [2,5,8,7,3]
pt1y = [4,8,6,9,9]
pt2y = [5,6,9,4,5]

t2 = absangle(pt1x,pt1y,pt2x,pt2y)
print(t2)
