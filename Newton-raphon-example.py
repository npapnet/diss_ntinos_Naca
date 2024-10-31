#%%

fx = lambda x: x**3 - x -1

fp = lambda x: 3*x**2 - 1

tol = 1e-666
exit_flag = False

x0 = 1.5
max_iter = 100
counter = 0
while not exit_flag:
    
    xnew = x0 - fx(x0)/fp(x0)
    print(f"x0 = {x0:10.7g}, f(x0) = {fx(x0):10.7g} , f'(x0) = {fp(x0):10.7g}, new x = {xnew:10.7g}")
    if abs(xnew-x0)<tol :
        exit_flag= True
        print('εχω σύγκλιση')
    
    counter += 1
    if counter>max_iter:
        print('Δεν έχω σύγκλιση')
        break
    x0 = xnew




