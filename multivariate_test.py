#!/usr/bin/env python3
# -*- coding: utf-8 -*-
''' to test mutivariate distribution generation'''
import matplotlib.pyplot as plt
import numpy as np


mean = [0, 0]
cov = [[1, 0], [0, 100]]

x, y = np.random.multivariate_normal(mean, cov, 5000).T
plt.plot(x, y, 'x')
plt.axis('equal')
plt.show()