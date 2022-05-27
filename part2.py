import matplotlib.pyplot as plt
import numpy as np

plt.dpi = 500

constants = {
    'GRAVITY': 8.65e-13,
    'MASS_OF_EARTH': 5.97e24,
    'MASS_OF_MOON': 7.35e22,
    'DISTANCE_MOON_TO_EARTH': 384400,
    'RADIUS_MOON': 1737.1,
    'LAGRANGE_POINT': 346007.85,
    'Gm_E': 5.16405000e12,
    'Gm_M': 6.3577500e10,
    'RADIUS_OF_EARTH': 6378}


class Function:
    def __init__(self, start_value, step_size, dydx, dydx_start, target_value):
        self.start_value = start_value
        self.step_size = step_size
        self.dydx = dydx
        self.dydx_start = dydx_start
        self.dydx_current_value = dydx_start
        self.target_value = target_value
        self.current_value = start_value

    def euler(self):
        self.current_value += self.step_size * self.dydx_current_value + 0.5 * (self.step_size ** 2) * self.dydx(
            self.current_value)

    # print(self.current_value)

    def velocity(self):
        self.dydx_current_value += self.step_size * self.dydx(self.current_value)

    def target_in_range(self):
        for i in range(1000000):
            if self.current_value >= self.target_value:
                return self.current_value, i / 1000, self.dydx_current_value, self.target_value, self.dydx(
                    self.current_value)
            self.euler()
            self.velocity()
        return -1


class VFunction:
    def __init__(self, start_valueX, start_valueY, step_size, d2ydx2_X, d2ydx2_Y, dydx_startX, dydx_startY,
                 secondOrder):
        # self.start_valueX = start_valueX
        # self.start_valueY = start_valueY
        self.start_value = np.array([start_valueX, start_valueY], dtype='float16')
        self.step_size = step_size
        self.d2ydx2_X = d2ydx2_X
        self.d2ydx2_Y = d2ydx2_Y
        self.dydx_startX = dydx_startX
        self.dydx_startY = dydx_startY
        self.dydx_current_valueX = dydx_startX
        self.dydx_current_valueY = dydx_startY
        self.current_valueX = self.start_value[0]
        self.current_valueY = self.start_value[1]
        self.points = np.array([[0.0, self.start_value[0], self.start_value[1]]])
        self.secondOrder = secondOrder

    def v_euler(self):
        second_order_x = 0
        second_order_y = 0
        if self.secondOrder > 0:
            second_order_x = 0.5 * self.d2ydx2_X(self.current_valueX, self.current_valueY) * self.step_size ** 2
            second_order_y = 0.5 * self.d2ydx2_Y(self.current_valueX, self.current_valueY) * self.step_size ** 2

        self.current_valueX += self.step_size * self.dydx_current_valueX + second_order_x
        self.current_valueY += self.step_size * self.dydx_current_valueY + second_order_y
        # print(self.current_value)

    def v_velocity(self):
        self.dydx_current_valueX += self.step_size * self.d2ydx2_X(self.current_valueX, self.current_valueY)
        self.dydx_current_valueY += self.step_size * self.d2ydx2_Y(self.current_valueX, self.current_valueY)

    def target_in_range(self):
        for i in range(500000):
            if self.current_valueX ** 2 + self.current_valueY ** 2 < 6370 ** 2:
                return self.current_valueX, i / 1000
            if i / 10 - np.modf(i / 10)[1] < 1 / 10:
                # self.points = np.insert(self.points,,[i/1000,self.current_valueX,self.current_valueY],axis=0)
                self.points = np.concatenate((self.points, [[i / 1000, self.current_valueX, self.current_valueY]]),
                                             axis=0)
            # print(self.points)
            # print(self.current_valueX)
            self.v_euler()
            self.v_velocity()
        return -1


def acceleration(x): return (constants['Gm_M'] / ((constants['DISTANCE_MOON_TO_EARTH'] - x) ** 2)
                             - constants['Gm_E'] / (x ** 2))


def acceleration_x(x, y): return constants['Gm_M'] * (constants['DISTANCE_MOON_TO_EARTH'] - x) / \
                                 pow((constants['DISTANCE_MOON_TO_EARTH'] - x) ** 2 + y ** 2, 3 / 2) \
                                 - constants['Gm_E'] * x / pow(x ** 2 + y ** 2, 3 / 2)


def acceleration_y(x, y): return -constants['Gm_M'] * y / pow((constants['DISTANCE_MOON_TO_EARTH'] - x) ** 2 + y ** 2,
                                                              3 / 2) \
                                 - constants['Gm_E'] * y / pow(x ** 2 + y ** 2, 3 / 2)


if __name__ == '__main__':
    GravityFunction = Function(6371, 0.001, acceleration, 40000,
                               (constants['DISTANCE_MOON_TO_EARTH'] - constants['RADIUS_MOON']))
    print(GravityFunction.target_in_range())
    SlingShot = VFunction(6563, 0, 0.001, acceleration_x, acceleration_y, 39210.0, -2100.0, 0)
    SlingShotSecondOrder = VFunction(6563, 0, 0.001, acceleration_x, acceleration_y, 39210.0, -2100.0, 1)
    # SlingShot2 = VFunction(6563,0,0.001,acceleration_x,acceleration_y,39280.0,-2310.0,0)
    # SlingShot3 = VFunction(6563,0,0.001,acceleration_x,acceleration_y,39215.0,-2210.0,0)
    # SlingShot4 = VFunction(6563.0,0,0.001,acceleration_x,acceleration_y,39297.0,-2097.0,0)
    # SlingShot5 = VFunction(6563.0,0,0.001,acceleration_x,acceleration_y,50000.0,-2000.0,0)
    print(SlingShot.target_in_range())
    print(SlingShotSecondOrder.target_in_range())
    # print(SlingShot2.target_in_range())
    # print(SlingShot3.target_in_range())
    # print(SlingShot4.target_in_range())
    # print(SlingShot5.target_in_range())
    float_formatter = "{:.3f}".format
    np.set_printoptions(formatter={'float_kind': float_formatter})
    # for x in range(len(SlingShot4.points)):
    #    print(SlingShot4.points[x])
    x = SlingShot.points[:, 1]
    y = SlingShot.points[:, 2]
    x2o = SlingShotSecondOrder.points[:, 1]
    y2o = SlingShotSecondOrder.points[:, 2]
    # x2 = SlingShot2.points[:,1]
    # y2 = SlingShot2.points[:,2]
    # x3 = SlingShot3.points[:,1]
    # y3 = SlingShot3.points[:,2]
    # x4 = SlingShot4.points[:,1]
    # y4 = SlingShot4.points[:,2]
    # x5 = SlingShot5.points[:,1]
    # y5 = SlingShot5.points[:,2]
    t = np.linspace(0, 6.3, 100)
    Ex = 6378 * np.cos(t)
    Ey = 6378 * np.sin(t)
    Mx = constants['DISTANCE_MOON_TO_EARTH'] + 1849 * np.cos(t)
    My = 1849 * np.sin(t)
    plt.plot(Mx, My)
    plt.plot(Ex, Ey)
    plt.plot(x, y)
    plt.plot(x2o, y2o)
    # plt.plot(x2,y2)
    # plt.plot(x3,y3)
    # plt.plot(x4,y4)
    # plt.plot(x5,y5)
    plt.ylim([-70000, 90000])
    plt.xlim([-42000, 500000])
    plt.title('Flight Path')
    plt.savefig('./Plots/FlightPath.png')
    plt.show()
