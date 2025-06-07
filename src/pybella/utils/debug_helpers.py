import matplotlib.pyplot as plt

def pl_sol(arr):
    plt.figure()
    plt.imshow(arr, origin="lower", aspect="auto")
    plt.colorbar()
    plt.show()