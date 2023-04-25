class Box:
    def __init__(self, xlo=0, xhi=0, ylo=0, yhi=0, zlo=0, zhi=0):
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

    @property
    def size(self):
        return (self.xhi - self.xlo, self.yhi - self.ylo, self.zhi - self.zlo)
