import matplotlib.pyplot as plt
plt.dpi = 500

constants = {
    'GRAVITY': 8.65e-13,
    'MASS_OF_EARTH': 5.97e24,
    'MASS_OF_MOON': 7.35e22,
    'DISTANCE_MOON_TO_EARTH': 384400,
    'RADIUS_MOON': 1737.1,
    'LAGRANGE_POINT': 346007.85
}


class Function:
    def __init__(self, start_value, step_size, dydx, dydx_start, target_value):
        self.start_value = start_value
        self.step_size = step_size
        self.dydx = dydx
        self.dydx_start = dydx_start
        self.dydx_current_value = dydx_start
        self.target_value = target_value
        self.current_value = start_value

        self.time_history = []
        self.displacement_history = []
        self.velocity_history = []
        self.acceleration_history = []

    def euler(self):
        self.displacement_history.append(self.current_value)
        self.current_value += self.step_size * self.dydx_current_value

    def velocity(self):
        self.velocity_history.append(self.dydx_current_value)
        self.acceleration_history.append(self.dydx(self.current_value))
        self.dydx_current_value += self.step_size * self.dydx(self.current_value)

    def target_in_range(self):
        for i in range(1000000):
            if self.current_value >= self.target_value:
                return self.current_value, i / 1000, self.dydx_current_value
            self.time_history.append(i / 1000)
            self.euler()
            self.velocity()
        return -1


def acceleration(x):
    return (constants['GRAVITY'] * (constants['MASS_OF_MOON'] / ((constants['DISTANCE_MOON_TO_EARTH'] - x) ** 2))) \
           - constants['GRAVITY'] * (constants['MASS_OF_EARTH'] / (x ** 2))


if __name__ == '__main__':
    GravityFunction = Function(6563, 0.001, acceleration, 39774.05, (constants['DISTANCE_MOON_TO_EARTH']-constants['RADIUS_MOON']))
    print(GravityFunction.target_in_range())
    print(GravityFunction.dydx_current_value)

    plt.figure(1)
    plt.plot(GravityFunction.time_history, GravityFunction.displacement_history)
    plt.xlabel('Time (h)')
    plt.ylabel('Displacement (km)')
    plt.title('Displacement')
    plt.savefig('./Plots/Displacement.png', bbox_inches='tight')

    plt.figure(2)
    plt.plot(GravityFunction.time_history, GravityFunction.velocity_history)
    plt.xlabel('Time (h)')
    plt.ylabel('Velocity (km/h)')
    plt.title('Velocity')
    plt.savefig('./Plots/Velocity.png', bbox_inches='tight')

    plt.figure(3)
    plt.plot(GravityFunction.time_history, GravityFunction.acceleration_history)
    plt.xlabel('Time (h)')
    plt.ylabel('Acceleration (km/h^2)')
    plt.title('Acceleration')
    plt.savefig('./Plots/Acceleration.png', bbox_inches='tight')