import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import os.path 
import pickle

def circular_stimuli (field = 250, frequancy = 4, blocked_field = 0, x0 = 0, y0 = 0, pattern = 'pos_expanding', centered = False):  
                    #        [um],          [hz],      [um] (radius),   [ms]
    
    time_cycle = 1 / float(frequancy) # [sec]

    # Circle should expand to cover the entire squared area. 
    # Ratio between bounding to bounded circles is sec(pi/4)=sqrt(2)
    effective_field = field * math.sqrt(2)

    expanding_radius  = lambda t: blocked_field + (t % time_cycle) * ((effective_field / 2 - blocked_field) / time_cycle)
    collapsing_radius = lambda t: ((effective_field / 2) - (t % time_cycle) * 
                                   ((effective_field / 2 - blocked_field) / time_cycle))
    
    # circle center point
    if centered:
        c_x0 = x0 
        c_y0 = y0
        x0 = c_x0 - (field / 2)
        y0 = c_x0 - (field / 2)
    else:
        c_x0 = x0 + (field / 2)
        c_y0 = y0 + (field / 2)
    
    #print(x0, y0, c_x0, c_y0)
    
    distance_point_circle = lambda x, y: math.sqrt((x - c_x0)**2 + (y - c_y0)**2)

    def is_inside_expanding (x, y, t):
        if within_bound(x, y):
            return blocked_field < distance_point_circle(x,y) < expanding_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_outside_expanding (x, y, t): 
        if within_bound(x, y):
            return distance_point_circle(x,y) > expanding_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_inside_collapsing (x, y, t):
        if within_bound(x, y):
            return blocked_field < distance_point_circle(x,y) < collapsing_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_outside_collapsing (x, y, t):
        if within_bound(x, y):
            return distance_point_circle(x,y) > collapsing_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def within_bound(x, y):
        if ((x0 + field > x > x0) and (y0 + field > y > y0)):
            return True

    if pattern == 'pos_expanding':
        return is_inside_expanding
    
    if pattern == 'neg_expanding':
        return is_outside_expanding
    
    if pattern == 'pos_collapsing':
        return is_inside_collapsing
    
    if pattern == 'neg_collapsing':
        return is_outside_collapsing

    return None
    
def alternating_expanding_circles (field = 250, frequancy = 4, blocked_field = 50, x0 = 0, y0 = 0, delay = 100, centered = True):  
                                 #        [um],          [hz],               [um] (radius)
    
    positive_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0,
                                       pattern = 'pos_expanding', centered = centered)
    
    negative_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'neg_expanding', centered = centered)
    
    time_cycle = (1 / float(frequancy)) * 1000 # [mSec]
    
    def is_activated (x, y, t):

        t = t - delay
        if t % (time_cycle * 2) < time_cycle:
            return positive_cycle (x, y, t % time_cycle)
        else:
            return negative_cycle (x, y, t % time_cycle)
        
    return is_activated

def alternating_collapsing_circles (field = 250, frequancy = 4, blocked_field = 50, x0 = 0, y0 = 0, delay = 100, centered = True):  
                                  #        [um],          [hz],               [um] (radius)
    
    positive_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'pos_collapsing', centered = centered)
    
    negative_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'neg_collapsing', centered = centered)
    
    time_cycle = (1 / float(frequancy)) * 1000 # [mSec]
    
    def is_activated (x, y, t):
        
        t = t - delay
        if t % (time_cycle * 2) < time_cycle:
            return negative_cycle (x, y, t % time_cycle)
        else:
            return positive_cycle (x, y, t % time_cycle)
        
    return is_activated

def drifting_bar (field_x = 500, bar_size_x = 100, velocity = 10, direction = 'R2L'): 
                                                   # pixels / sec
    def is_activated (x, y, t):
        
        if direction == 'R2L':
            x_loc = t * velocity
            x_s_loc = x_loc - bar_size_x
        else:
            x_s_loc = field_x - (t * velocity)
            x_loc = x_s_loc + bar_size_x
              
        if x_loc > x > x_s_loc:
            return True
        else:
            return False

    return is_activated   
  
def alternating_bar (field_x = 500, bar_size_x = 100, velocity = 10):
    
    t_to_finish = (field_x + bar_size_x) / velocity
    
    left_moving_bar  = drifting_bar(field_x = field_x, bar_size_x = bar_size_x, velocity = velocity, direction = 'R2L')
    right_moving_bar = drifting_bar(field_x = field_x, bar_size_x = bar_size_x, velocity = velocity, direction = 'L2R')
    
    def is_activated (x, y, t):
        
        t = t % (t_to_finish * 2)
        
        if t < t_to_finish:
            return left_moving_bar (x, y, t)
        else:
            return right_moving_bar (x, y, t % t_to_finish)
        
    return is_activated
    
        
# *******************************************
# ********** Evaluation method  *************
# *******************************************

def evaluate_stimuli_pattern (t,      stimuli,            field_x, field_y): 
                           # [mSec],  [function pointer], [um] 

    res = np.zeros((field_y, field_x))
    
    for x in range(field_x):
        for y in range(field_y):
            if stimuli(x, y, t):
                res[y][x] = 1
    return res


    


