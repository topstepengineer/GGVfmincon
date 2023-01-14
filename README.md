# GGVfmincon
matlab code to generate GGV from a bicycle model. fmincon is used to find max accx for given accy.
mybike.m first finds the max lateral acceleration (acc_y) with zero longitudinal acceleration (acc_x) at a given speed. Then tries to find max acc_x for a range of acc_y.
