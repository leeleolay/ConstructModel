def compute_scaling_factors(D0_M, D0_P, Y_M, Y_P, eta0_M, eta0_P):
    # 计算长度的scaling factor
    lambda_ = D0_P / D0_M  # meters

    # 计算能量的scaling factor
    xi = (Y_P / Y_M) * (D0_P / D0_M)**2  # Joules

    # 计算力的scaling factor
    epsilon = (Y_P / Y_M) * (D0_P / D0_M)  # Newton

    # 计算压力的scaling factor
    nu = epsilon / lambda_**3  # Pascal

    # 计算时间的scaling factor
    tau = (eta0_P / eta0_M) * (Y_M / Y_P) * (D0_P / D0_M)  # seconds

    # 计算质量的scaling factor
    mu = (Y_P / Y_M) * tau**2  # kg

    return lambda_, xi, epsilon, nu, tau, mu

# 测试计算的函数
D0_P = 7.82E-6  # Meters
Y_P = 18.9E-6   # Newton/Meter = Joules/(Meter^2) for 2-Dimensions
eta0_P = 0.022  # Pascal*Second

D0_M =  8.25 # lambda
Y_M = 400   # epsilon/lambda=xi/(lambda^2)  for 2-Dimensions
eta0_M = 168.87  # nu*tau

# 使用函数进行计算
lambda_, xi, epsilon, nu, tau, mu = compute_scaling_factors(D0_M, D0_P, Y_M, Y_P, eta0_M, eta0_P)

print("Length scaling factor (lambda):", '{:.8e}'.format(lambda_))
print("Energy scaling factor (xi):", '{:.8e}'.format(xi))
print("Force scaling factor (epsilon):", '{:.8e}'.format(epsilon))
print("Pressure scaling factor (nu):", '{:.8e}'.format(nu))
print("Time scaling factor (tau):", '{:.8e}'.format(tau))
print("Mass scaling factor (mu):", '{:.8e}'.format(mu))