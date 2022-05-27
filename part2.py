import math
import matplotlib.pyplot as plt

from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

plt.dpi = 500

G = 8.65e-13
M_EARTH = 5.97e24
M_MOON = 7.35e22
D_ME = 384400
R_MOON = 1737.1
LAGRANGE_POINT = 346007.85


class Function:
    def __init__(self, start_value_x, start_value_y, step_size, dydx_x, dydx_y, dydx_start_x, dydx_start_y,
                 target_value_x, target_value_y):
        self.start_value_x = start_value_x
        self.start_value_y = start_value_y
        self.step_size = step_size
        self.dydx_x = dydx_x
        self.dydx_y = dydx_y
        self.dydx_start_x = dydx_start_x
        self.dydx_start_y = dydx_start_y
        self.target_value_x = target_value_x
        self.target_value_y = target_value_y

        self.current_value_x = self.start_value_x
        self.current_value_y = self.start_value_y

        self.dydx_current_value_x = self.dydx_start_x
        self.dydx_current_value_y = self.dydx_start_y

        self.time_history = []
        self.displacement_history_x = []
        self.displacement_history_y = []
        self.displacement_history_total = []
        self.velocity_history = []
        self.acceleration_history = []
        self.max_displacement_x = 0

    def euler(self):
        self.displacement_history_total.append(math.sqrt(self.current_value_x ** 2 + self.current_value_y ** 2))
        self.current_value_x += self.step_size * self.dydx_current_value_x
        self.current_value_y += self.step_size * self.dydx_current_value_y

    def velocity(self):
        self.velocity_history.append(math.sqrt(self.dydx_current_value_x ** 2 + self.dydx_current_value_y ** 2))
        self.acceleration_history.append(math.sqrt(
            self.dydx_x(self.current_value_x, self.current_value_y) ** 2 + self.dydx_y(self.current_value_x,
                                                                                       self.current_value_y) ** 2))
        self.dydx_current_value_x += self.step_size * self.dydx_x(self.current_value_x, self.current_value_y)
        self.dydx_current_value_y += self.step_size * self.dydx_y(self.current_value_x, self.current_value_y)

    def target_in_range(self):
        for i in range(250000):
            self.displacement_history_x.append(self.current_value_x)
            self.displacement_history_y.append(self.current_value_y)
            if self.current_value_x > self.max_displacement_x:
                self.max_displacement_x = self.current_value_x
            if self.current_value_x <= self.target_value_x and self.current_value_y <= self.target_value_y and self.max_displacement_x >= D_ME:  # TODO: updated displacement history. Make it into two lists
                return self.current_value_x, self.current_value_y, i / 1000, self.dydx_current_value_x, self.dydx_current_value_y
            self.time_history.append(i / 1000)
            self.euler()
            self.velocity()
        # return -1


def acceleration_x(x, y):
    return ((-G * M_EARTH * x) / (
            math.sqrt((x ** 2) + (y ** 2)) * ((x ** 2) + (y ** 2)))) + (((D_ME - x) *
                                                                         G * M_MOON) / (math.sqrt(
        (D_ME - x) ** 2 + (y ** 2)) * (((D_ME - x) ** 2) + (y ** 2))))


def acceleration_y(x, y):
    return ((-G * M_EARTH * y) / (
            math.sqrt(x ** 2 + y ** 2) * (x ** 2 + y ** 2))) + ((y *
                                                                 -G * M_MOON) / (
                                                                        math.sqrt((D_ME - x) ** 2 + (y ** 2)) * (
                                                                        ((D_ME - x) ** 2) + (y ** 2))))


if __name__ == '__main__':
    f = Function(6563, 0, 0.001, acceleration_x, acceleration_y, 39230, -1275, 6563, 0)
    print(f.target_in_range())

    plt.figure(1)
    plt.plot(f.displacement_history_x, f.displacement_history_y)
    circle1 = plt.Circle((0, 0), 6563, color='orange', fill=False)
    circle2 = plt.Circle((D_ME, 0), R_MOON, color='r', fill=False)
    plt.xlim(-20000, 500000)
    plt.ylim(-65000, 80000)
    plt.ylabel('Displacement (km)')
    plt.ylabel('Displacement (km)')
    plt.title('Flight Path')

    plt.annotate('1', xy=(204996.8327792941, -11495.180125773899))
    plt.annotate('2', xy=(292228.0650928426, -17148.717065768455))
    plt.annotate('3', xy=(347977.8642659397, -19444.202147591983))
    plt.annotate('4', xy=(380949.7141145375, 11086.064512151816))
    plt.annotate('5', xy=(328582.2212321744, 12989.859997015128))
    plt.annotate('6', xy=(263400.4594528887, 6491.854550348502))
    plt.annotate('7', xy=(154465.5669197549, -1233.7123160741326))
    plt.annotate('8', xy=(6546.932540818372, -2545.504565159628))

    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    plt.savefig('./Plots/FlightPath.png')
    plt.show()

    plt.figure(2)
    plt.plot(f.time_history, f.velocity_history)
    plt.xlabel('Time (h)')
    plt.ylabel('Velocity (km/h)')
    plt.title('Velocity2D')
    plt.savefig('./Plots/Velocity2D.png', bbox_inches='tight')
    plt.show()

    plt.figure(3)
    plt.plot(f.time_history, f.acceleration_history)
    plt.xlabel('Time (h)')
    plt.ylabel('Acceleration (km/h^2)')
    plt.title('Acceleration2D')
    plt.savefig('./Plots/Acceleration2D.png', bbox_inches='tight')
    plt.show()

    plt.figure(4)
    plt.plot(f.time_history, f.displacement_history_total)
    plt.xlabel('Time (h)')
    plt.ylabel('Displacement (km)')
    plt.title('Displacement2D')
    plt.savefig('./Plots/Displacement2D.png', bbox_inches='tight')
