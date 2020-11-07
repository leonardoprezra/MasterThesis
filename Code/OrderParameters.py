import freud
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# Try to plot using KDE if available, otherwise revert to histogram
try:
    from sklearn.neighbors.kde import KernelDensity
    kde = True
except:
    kde = False

np.random.seed(1)
