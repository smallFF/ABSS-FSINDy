from pprint import pprint

def generate_lorenz96_model(n, F):
    # file_name = 'lorenz96_m_'+str(m)+'_F_'+str(F)
    def ode(t, x):
        dx = []
        # formated_str = []
        for i in range(n):
            ind_i_plus_1 = modify_ind(i+1, n)
            ind_i_minus_2 = modify_ind(i-2, n)
            ind_i_minus_1 = modify_ind(i-1, n)
            dx.append(
                (x[ind_i_plus_1] - x[ind_i_minus_2])*x[ind_i_minus_1] - x[i] + F       
            )
        #     formated_str.append(
        #         f'dX[{i}] = (X[{ind_i_plus_1}] - X[{ind_i_minus_2}])*X[{ind_i_minus_1}] - X[{i}] + {F}'
        #     )
        # pprint(formated_str)
        return dx
    return ode
        
        
def modify_ind(i, n):
    if i >=n:
        ind = i - n
    elif i < 0:
        ind = i + n
    else:
        ind = i
    
    return ind;

if __name__ == '__main__':
    from generate_lorenz96_model import generate_lorenz96_model
    n = 6
    F = 8
    rhs = generate_lorenz96_model(n, F)
    x0 = [8 for i in range(n)]
    x0[0] = 1
    pprint(rhs(0, x0))