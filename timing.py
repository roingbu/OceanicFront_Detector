
import time
from time import clock as clk

def print_timing(func):
    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        #print "%s took %f ms with args: %s" % (func.func_name, (t2-t1)*1000.0, arg)
        print "%s took %f ms" % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper

#----------------------------------------------------------------------

def timing_collector(out):
    while True:
        (name, t) = (yield)
        if name in out:
            out[name] += [t]
        else:
            out[name] = [t]


class Timing(object):
    def __init__(self):#,cor):
        self.timings = {}
        self.col = self.__collector()
        self.col.next()

    def __collector(self):
        while True:
            (name, t) = (yield)
            if name in self.timings:
                self.timings[name] += [t]
            else:
                self.timings[name] = [t]

    def __call__(self, func):
        def wrapper(*arg, **kwargs):
            t1 = time.time()
            res = func(*arg, **kwargs)
            t2 = time.time()
            t = (t2-t1)*1000.0  #time in milliseconds
            data = (func.__name__, t)
            self.col.send(data)
            return res
        return wrapper

    def __str__(self):
        return "%s" % self.timings

#----------------------------------------------------------------------


class Timer:
    def __init__(self, msg):
        self.msg = msg
        self.running_time = 0.

    def tic(self):
        self.time = clk()

    def toc(self, out = 0):
        self.time = clk() - self.time
        self.running_time += self.time
        if out != 0:
            if self.msg:
                print "(%s), time: %f sec" % (self.msg, self.time)
            else:
                print "time: %f sec" % (self.time)

    def reset(self):
        self.running_time = 0.

    def xprint(self):
        if self.msg:
            print "(%s), time: %f sec" % (self.msg, self.running_time)
        else:
            print "time: %f sec" % (self.running_time)


#----------------------------------------------------------------------

if __name__ == "__main__":

    @print_timing
    def add(x,y):
        for i in range(10000):
            c = x + y
        return c


    add(3.,4.)

    t = Timer()
    x = 3.

    t.tic()
    add(5.,7.)
    t.toc()

    t.tic()
    for i in range(1000):
        #x = log(2.+cos(x))
        x = 3
    t.toc()


