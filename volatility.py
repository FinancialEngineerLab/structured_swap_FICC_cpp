import QuantLib as ql
import numpy as np
import pandas as pd



# settlementDate = ql.Date(31, ql.December, 2021)
settlementDate = ql.Date(30, ql.April, 2021)
calendar = ql.SouthKorea()



class StructuredSwap:
    def __init__(self, maturityDate, indexes) -> None:
        self.maturityDate = maturityDate
        self.indexes = indexes


def period_to_years(p: ql.Period):
    if p.units() == ql.Months:
        return p.length() / 12.0
    elif p.units() == ql.Years:
        return p.length()


def index_curve():
    zeroRates = [('0D',0.00489990232726306),
                ('1M',0.00489990232726306),
                ('3M',0.00729336504901711),
                ('6M',0.0072183976819218),
                ('9M',0.00756841947715215),
                ('1Y',0.00811894680946977),
                ('18M',0.0094219812222739),
                ('2Y',0.010451343656196),
                ('3Y',0.0124436447952613),
                ('4Y',0.0140421854966705),
                ('5Y',0.0150603799469474),
                ('6Y',0.0157229797860345),
                ('7Y',0.0162608367063139),
                ('8Y',0.016618275595913),
                ('9Y',0.017087893685427),
                ('10Y',0.0173991409047416),
                ('12Y',0.0178673827930665),
                ('15Y',0.0176907924545021),
                ('20Y',0.0175718994031361),
                ('25Y',0.0177055936587631),
                ('30Y',0.0176898158517235)]

    dates = []
    rates = []

    for rate in zeroRates:
        d = calendar.advance(settlementDate, ql.Period(rate[0]), ql.ModifiedFollowing)
        dates.append(d)
        rates.append(rate[1])

    curve = ql.ZeroCurve(dates, rates, ql.Actual365Fixed(), calendar)

    return ql.YieldTermStructureHandle(curve)


def vol_matrix():
    atmOptionTenors = [
        ql.Period(1, ql.Months),
        ql.Period(3, ql.Months),
        ql.Period(6, ql.Months),
        ql.Period(9, ql.Months),
        ql.Period(1, ql.Years),
        ql.Period(2, ql.Years),
        ql.Period(3, ql.Years),
        ql.Period(18, ql.Months),
        ql.Period(10, ql.Years),
    ]

    atmSwapTenors = [ql.Period(y, ql.Years) for y in [1,2,3,4,5,7,10]]

    m = np.array(
            [[3200, 3963, 4725, 5488, 6250, 6225, 6188],
             [4000, 4356, 4712, 5069, 5425, 5628, 5933],
             [4041, 4341, 4641, 4942, 5242, 5466, 5802],
             [4082, 4313, 4545, 4776, 5008, 5273, 5670],
             [4125, 4310, 4495, 4680, 4865, 5132, 5532],
             [4414, 4513, 4613, 4713, 4813, 5020, 5331],
             [4727, 4734, 4742, 4749, 4756, 4899, 5112],
             [4339, 4386, 4433, 4480, 4528, 4571, 4635],
        ]
    ) / 1000000

    atmVol = ql.SwaptionVolatilityStructureHandle(
        ql.SwaptionVolatilityMatrix(
            calendar,
            ql.Following,
            atmOptionTenors,
            atmSwapTenors,
            ql.Matrix(m.tolist()),
            ql.Actual365Fixed(),
        )
    )

    atmVol.atmOptionTenors = atmOptionTenors
    atmVol.atmSwapTenors = atmSwapTenors
    atmVol.m = m

    return atmVol


def vol_matrix_file():
    df = pd.read_excel('Market.xlsx', sheet_name=2, index_col=0)

    atmOptionTenors = [ql.Period(p) for p in df.index]
    atmSwapTenors = [ql.Period(p) for p in df.columns]

    m = df.to_numpy()

    atmVol = ql.SwaptionVolatilityStructureHandle(
        ql.SwaptionVolatilityMatrix(
            calendar,
            ql.ModifiedFollowing,
            atmOptionTenors,
            atmSwapTenors,
            ql.Matrix(m.tolist()),
            ql.Actual365Fixed(),
        )
    )

    atmVol.atmOptionTenors = atmOptionTenors
    atmVol.atmSwapTenors = atmSwapTenors
    atmVol.m = m
    atmVol.df = df
    
    x = [ period_to_years(tenor) for tenor in atmSwapTenors ]
    y = [ period_to_years(tenor) for tenor in atmOptionTenors ]
    
    print(x)
    print(y)

    atmVol.interpolation = ql.BilinearInterpolation(x, y, m.tolist())

    return atmVol


def swaption_helper(exerciseTenor, swapTenor, vol):
    fixedLegTenor = ql.Period(3, ql.Months)
    dayCounter = ql.Actual365Fixed()
    curve = index_curve()
    index = ql.Libor('libor', ql.Period(3, ql.Months), 1, ql.KRWCurrency(), ql.SouthKorea(), dayCounter, curve)
    volatilityType = ql.Normal
    
    helper = ql.SwaptionHelper(exerciseTenor, 
                                swapTenor,
                                ql.QuoteHandle(ql.SimpleQuote(vol)),
                                index,
                                fixedLegTenor,
                                dayCounter,
                                dayCounter,
                                curve,
                                ql.BlackCalibrationHelper.RelativePriceError,
                                ql.nullDouble(),
                                1.0,
                                volatilityType)
    helper.vol = vol

    return helper


def swaptions_interpolate(volMatrix, x, exerciseDate, y, endDate, curve):
    dayCounter = ql.Actual365Fixed()

    # x = dayCounter.yearFraction(exerciseDate, endDate)
    # y = dayCounter.yearFraction(settlementDate, exerciseDate)
    vol = volMatrix.interpolation(y, x)
    
    fixedLegTenor = ql.Period(3, ql.Months)
    index = ql.Libor('libor', ql.Period(3, ql.Months), 1, ql.KRWCurrency(), ql.SouthKorea(), dayCounter, curve)
    volatilityType = ql.Normal

    helper = ql.SwaptionHelper(exerciseDate,
                                endDate,
                                ql.QuoteHandle(ql.SimpleQuote(vol)),
                                index,
                                fixedLegTenor,
                                dayCounter,
                                dayCounter,
                                curve,
                                ql.BlackCalibrationHelper.RelativePriceError,
                                ql.nullDouble(),
                                1.0,
                                volatilityType)
                                
    helper.vol = vol

    return helper


def swaptions_diagonal(swap: StructuredSwap):
    dayCounter = ql.Actual365Fixed()
    volMatrix = vol_matrix_file()
    maturity_t = dayCounter.yearFraction(settlementDate, swap.maturityDate)
    curve = index_curve()

    swap_t = min(int(maturity_t), 10)
    option_t = maturity_t - swap_t

    instruments = []
    print(option_t, maturity_t)
    for i in range(0, 10):
        roop_option_t = option_t + i
        roop_swap_t = swap_t - i

        if roop_option_t <= 10 and swap_t - i > 0:
            pair = (roop_option_t, roop_swap_t)
            exerciseDate = ql.NullCalendar().advance(settlementDate, ql.Period(int(roop_option_t*365), ql.Days) )
            swapMaturityDate = calendar.advance(exerciseDate, ql.Period(roop_swap_t, ql.Years))
            
            helper = swaptions_interpolate(volMatrix, pair[0], exerciseDate, pair[1], swapMaturityDate, curve)
            # print(pair, ql.Actual365Fixed().yearFraction(settlementDate, exerciseDate), ql.Actual365Fixed().yearFraction(settlementDate, swapMaturityDate), exerciseDate.ISO(), swapMaturityDate.ISO(), helper.vol)

            print(pair, helper.vol)
            instruments.append(helper)

    print(len(instruments))
    return instruments


def swaptions_diagonal_and_index(swap: StructuredSwap):
    dayCounter = ql.Actual365Fixed()
    volMatrix = vol_matrix_file()
    maturity_t = dayCounter.yearFraction(settlementDate, swap.maturityDate)
    curve = index_curve()

    swap_t = min(int(maturity_t), 10)
    option_t = maturity_t - swap_t

    instruments = []
    print(option_t, maturity_t)
    for i in range(0, 10):
        roop_option_t = option_t + i
        roop_swap_t = swap_t - i

        if roop_option_t <= 10 and swap_t - i > 0:
            pair = (roop_option_t, roop_swap_t)
            exerciseDate = ql.NullCalendar().advance(settlementDate, ql.Period(int(roop_option_t*365), ql.Days) )
            swapMaturityDate = calendar.advance(exerciseDate, ql.Period(roop_swap_t, ql.Years))
            
            helper = swaptions_interpolate(volMatrix, pair[0], exerciseDate, pair[1], swapMaturityDate, curve)
            # print(pair, ql.Actual365Fixed().yearFraction(settlementDate, exerciseDate), ql.Actual365Fixed().yearFraction(settlementDate, swapMaturityDate), exerciseDate.ISO(), swapMaturityDate.ISO(), helper.vol)
            print(pair[0], pair[1], helper.vol)
            # print(pair, helper.swaptionExpiryDate(), helper.vol)
            
            instruments.append(helper)
    
    print(ql.Settings.instance().evaluationDate)

    for swap_maturity in swap.indexes:
        for exerciseTenor in volMatrix.atmOptionTenors:
            # print(type(exerciseTenor), type(swap_maturity))
            vol = volMatrix.df[str(swap_maturity).lower()][str(exerciseTenor).lower()]
            helper = swaption_helper(exerciseTenor, swap_maturity, vol)

            # print(exerciseTenor, swap_maturity, helper.swaptionExpiryDate(), helper.vol)
            # print(exerciseTenor, swap_maturity, helper.vol)
            print(period_to_years(exerciseTenor), period_to_years(swap_maturity), helper.vol)
            
            instruments.append(helper)
    

    print(len(instruments))
    return instruments



def swaptions_fullset():
    atmVol = vol_matrix_file()
    curve = index_curve()

    fixedLegTenor = ql.Period(3, ql.Months)
    dayCounter = ql.Actual365Fixed()
    index = ql.Libor('libor', ql.Period(3, ql.Months), 1, ql.KRWCurrency(), ql.SouthKorea(), dayCounter, curve)
    volatilityType = ql.Normal

    helpers = []

    for r, exerciseTenor in enumerate(atmVol.atmOptionTenors):
        for c, swapTenor in enumerate(atmVol.atmSwapTenors):
            helper = ql.SwaptionHelper(exerciseTenor, 
                                    swapTenor,
                                    ql.QuoteHandle(ql.SimpleQuote(atmVol.m[r][c])),
                                    index,
                                    fixedLegTenor,
                                    dayCounter,
                                    dayCounter,
                                    curve,
                                    ql.BlackCalibrationHelper.RelativePriceError,
                                    ql.nullDouble(),
                                    1.0,
                                    volatilityType)

            helpers.append(helper)


    return helpers

    # normal to 


class ModelCalibrator:
    def __init__(self, endCriteria):        
        self.endCriteria = endCriteria
        self.helpers = []

    def AddCalibrationHelper(self, helper):
        self.helpers.append(helper)

    def Calibrate(self, model, engine, curve, fixedParameters):
        # assign pricing engine to all calibration helpers
        for i in range(len(self.helpers)):
            self.helpers[i].setPricingEngine(engine)
        method = ql.LevenbergMarquardt()
        if(len(fixedParameters) == 0):
            print(len(self.helpers))
            print(method)
            print(self.endCriteria)
            # model.calibrate(self.helpers, method, self.endCriteria)
            model.calibrate(self.helpers, method, self.endCriteria)
        else:
            model.calibrate(self.helpers, method, self.endCriteria,
                ql.NoConstraint(), [], fixedParameters)


def main():

    # general parameters
    ql.Settings.instance().evaluationDate = settlementDate

    # create market data: term structure and diagonal volatilities
    curveHandle = index_curve()
    # vol = ql.CreateSwaptionVolatilityList()

    # create calibrator object
    endCriteria = ql.EndCriteria(10000, 100, 0.000001, 0.00000001, 0.00000001)
    calibrator = ModelCalibrator(endCriteria)
    
    # SWP16060016014
    maturityDate = ql.DateParser_parseISO('2031-01-17') # 2031-02-17
    indexes = [ql.Period('2y'), ql.Period('10y')]
    swap = StructuredSwap(maturityDate, indexes)

    # for h in swaptions_fullset():
    # for h in swaptions_diagonal(swap):
    for h in swaptions_diagonal_and_index(swap):
        calibrator.AddCalibrationHelper(h)

    # create model and pricing engine, calibrate model and print calibrated parameters
    model = ql.G2(curveHandle, a=0.05, b=0.5, rho=-0.9)
    engine = ql.G2SwaptionEngine(model, range=2, intervals=100)

    fixedParameters = [True, False, True, False, True]
    # fixedParameters = [False, False, False, False, False]
    calibrator.Calibrate(model, engine, curveHandle, fixedParameters)
    # a=0.1, sigma=0.01, b=0.1, eta=0.01, rho=-0.75
    print(model.params())


main()

# v = vol_matrix_file()
# print(v.df['2y']['3m'])
# print(v.interpolation(2,5))

