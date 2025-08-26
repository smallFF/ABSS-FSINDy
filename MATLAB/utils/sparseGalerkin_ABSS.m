function dx = sparseGalerkin_ABSS(t,X, max_polyorder, Xi)
dx = (getLibrary(X', max_polyorder)*Xi)';