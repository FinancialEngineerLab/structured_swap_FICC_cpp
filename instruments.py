import QuantLib as ql

# 자... 지금 우선 하려고 하는 것은.
# 구조화 스왑의 평가 를 여기서...?! 음... 머.. 속도는... 어... 느리면 머...
# 가능? pricing engine 을

# 지금 이게 왜 이런 마음이 자꾸 드냐면,,, mxdevtool 이거가 음... 자꾸 마음에 걸려. 자꾸 완벽하게 하려고 하는 거 같고,
# 그래서리. 이게 쌓이지 않는다는거.
# 쌓이려면 소스가 공개되어야함. -> 이 상태로 으럅 gen 하는 것. 시나리오 여기서는 하나로 두개를 뽑아내는 거임
# 이 코드를 가지고 나중에 거기다가 마이그레이션 시켜야할 듯. 음...
# 이게 엔진이니까... 음... 자꾸 파이썬으로 되어서리.
# 머가 효율적인지 모르겠네
# 음... 글게.. 이런.. 이게 무슨.. 허허허.. 아 그거 돈좀 벌었음.


class StructuredLeg:
    def __init__(self) -> None:
        pass


class CMSSpreadRange:
    def __init__(self) -> None:
        pass


class StrucutredSwap(ql.Instrument):
    def __init__(self, payLeg=None, recLeg=None, exerciseLeg=None) -> None:
        self.payLeg = payLeg
        self.recLeg = recLeg
        self.exerciseLeg = exerciseLeg
        self.maturityDate = ql.Date(30, ql.April, 2032)
        self.engine = None
        self.results = dict()

    def setPricingEngine(self, arg2):
        print('setPricingEngine')
        self.engine = arg2

    def NPV(self):
        self.engine.calculate(self)
        return self.results['NPV']

    def greeks(self):
        return self.results['greeks']



class MCStrucutredSwapEngine(ql.PricingEngine):
    def __init__(self, process, timeStep, samples, seed, evaluationDate) -> None:
        self.process = process
        self.timeStep = timeStep
        self.samples = samples
        self.seed = seed
        self.evaluationDate = evaluationDate


    def calculate(self, instrument):
        fixedLeg = instrument.payLeg
        floatingLeg = instrument.recLeg

        maturity_t = ql.Actual365Fixed().yearFraction(self.evaluationDate, instrument.maturityDate)
        timegrid = ql.TimeGrid(maturity_t, self.timeStep)

        rng = ql.GaussianRandomSequenceGenerator(
            ql.UniformRandomSequenceGenerator(self.timeStep, ql.UniformRandomGenerator(seed=1)))

        # seq = ql.GaussianPathGenerator(self.process, 10, self.timeStep, rng, False)
        seq = ql.GaussianPathGenerator(self.process, timegrid, rng, False)
        path = seq.next().value()
        print(path[1])

        # 이게 이런식으로 되면 나도 하기 쉽고, 추가하는 거도 쉽고,
        # 추가가 필요하면, 으럇.

        instrument.results['NPV'] = 0.1





evaluationDate = ql.Date(30, ql.April, 2022)

swap = StrucutredSwap(None, None, None)

spot_curve_handle = ql.YieldTermStructureHandle(ql.FlatForward(evaluationDate, 0.02, ql.Actual365Fixed()))
hw_process = ql.HullWhiteProcess(spot_curve_handle, a=0.05, sigma=0.01)
engine = MCStrucutredSwapEngine(hw_process, 365, 1000, 1, evaluationDate)


swap.setPricingEngine(engine)
npv = swap.NPV()
print(npv)

