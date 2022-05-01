#pragma once

#include <ql/instruments/swap.hpp>
#include <ql/exercise.hpp>
#include <ql/pricingengines/mclongstaffschwartzengine.hpp>
#include <ql/cashflows/all.hpp>
#include <ql/processes/hullwhiteprocess.hpp>
#include <ql/processes/hullwhite1factorprocess.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/methods/montecarlo/lsmbasissystem.hpp>
#include <ql/indexes/ibor/all.hpp>
#include <ql/time/calendars/all.hpp>
#include <qle/termstructures/yieldcurveext.hpp>
#include <ql/cashflows/couponpricer.hpp>
#include <ql/time/daycounters/all.hpp>
#include <ql/models/parameter.hpp>
#include <ql/math/interpolations/all.hpp>

using namespace QuantLib;

namespace QuantLib
{
	class LegExerciseOption {
	public:
		LegExerciseOption(const std::vector<Date>& dates,
			const std::vector<Date>& settlementDates)
			: dates_(dates), settlementDates_(settlementDates)
		{
		}

		const std::vector<Date>& dates() const { return dates_; }
		Date lastDate() const { return dates_.back(); }

	protected:
		std::vector<Date> dates_;
		std::vector<Date> settlementDates_;
	};


	class LinearInterpolationParameter : public Parameter {
	private:
		class Impl : public Parameter::Impl {
		public:
			Impl(const std::vector<Time>& times)
				: times_(times) {}

			Real value(const Array& params, Time t) const {

				if (t <= times_[0])
					return params[0];

				if (times_.back() < t)
					return params.back();

				Size size = times_.size();
				for (Size i = 1; i < size; i++) {

					if (t <= times_[i])
					{
						Real a = (params[i] - params[i - 1]);
						Real b = (times_[i] - times_[i - 1]);
						Real c = (t - times_[i - 1]);

						return params[i-1] + (t - times_[i - 1]) * (params[i] - params[i-1]) / (times_[i] - times_[i - 1]);
					}

				}
				return params[size-1];
			}
			std::vector<Time> times() { return times_; }


		private:
			std::vector<Time> times_;
		};
	public:
		LinearInterpolationParameter(const std::vector<Time>& times,
			const std::vector<Real>& values,
			const Constraint& constraint = NoConstraint())
			: Parameter(times.size(),
				boost::shared_ptr<Parameter::Impl>(
					new LinearInterpolationParameter::Impl(times)), constraint)
		{
			for (Size i = 0; i < times.size(); i++)
				this->setParam(i, values[i]);

		}
	};
}

namespace TestEngine
{
	class ZeroCouponCallableSwap : public Swap
	{
	public:
		class arguments : public virtual PricingEngine::arguments {
		public:
			Leg fixedLeg;
			Leg floatingLeg;
			std::vector<Date> excerciseDates;
			void validate() const {}
		};

		class results : public Instrument::results {
		public:
			void reset() {}
		};

		class engine : public GenericEngine<ZeroCouponCallableSwap::arguments,
			ZeroCouponCallableSwap::results> {};
		//ZeroCouponCallableSwap(
		//	std::vector<Leg> legs,
		//	const std::vector<bool>& payer,
		//	const std::vector<boost::shared_ptr<LegExerciseOption>> options)
		//: Swap(legs, payer), options_(options)
		//{
		//}

		// zero callable 이니까 fixed 기준으로 함.
		// options 에 들어있는 ex date는 floatingLeg의 subset이어야함?
		ZeroCouponCallableSwap(
			Leg fixedLeg,
			Leg floatingLeg,
			const boost::shared_ptr<LegExerciseOption>& option,
			bool isFixedRec)
			: Swap(fixedLeg, floatingLeg), option_(option)
		{
			boost::shared_ptr<IborCouponPricer> pricer(new BlackIborCouponPricer);
			setCouponPricer(floatingLeg, pricer);
		}

		void setupArguments(PricingEngine::arguments* args) const
		{
			ZeroCouponCallableSwap::arguments* arguments = dynamic_cast<ZeroCouponCallableSwap::arguments*>(args);
			QL_REQUIRE(arguments != 0, "wrong argument type");

			arguments->excerciseDates = this->option_->dates();
			arguments->fixedLeg = this->leg(0);
			arguments->floatingLeg= this->leg(1);
		}

		void fetchResults(const PricingEngine::results*) const
		{

		}

	protected:
		boost::shared_ptr<LegExerciseOption> option_;

	};

}


namespace QuantLib
{
	// -------------------------------------------------------------------------------------------

	//template <class RNG = PseudoRandom, class S = Statistics>
	//class MCZeroCouponCallableSwapEngine
	//	: public MCLongstaffSchwartzEngine<TestEngine::ZeroCouponCallableSwap::engine,
	//	SingleVariate, RNG, S> {
	class MCZeroCouponCallableSwapEngine2 : public TestEngine::ZeroCouponCallableSwap::engine {
	  public:
		MCZeroCouponCallableSwapEngine2(
			//const boost::shared_ptr<StochasticProcess>& process,
			const boost::shared_ptr<HullWhite1FactorProcess>& process,
			Size timeSteps,
			Size samples,
			BigNatural seed,
			Size nCalibrationSamples)
		: process_(process),
		  timeSteps_(timeSteps), samples_(samples), seed_(seed), nCalibrationSamples_(nCalibrationSamples)
		{
			this->hw_model_ = boost::shared_ptr<HullWhite>(
				new HullWhite(process->fitting(), process->a(10.0), process->sigma(10.0)));
		}

		Real legNpv(const Leg& leg, const Path& path, Time t)
		{
			return 0.0;
		}

		Real rate(const Path& path, Time t) const
		{
			Size pos = path.timeGrid().closestIndex(t, -1);

			Real r1 = path[pos];
			Real r2 = path[pos+1];
			Real t1 = path.timeGrid().at(pos);
			Real t2 = path.timeGrid().at(pos+1);
			Real res = r1 * (t2 - t) + r2 * (t - t1);

			return res / (t2 - t1);

		}

		//Real compound(const Path& path, Date from, Date to, DayCounter dc) const
		//{
		//	Date evaluationDate = Settings::instance().evaluationDate();

		//	Real fraction = dc.yearFraction(from, to);
		//	Real from_t = dc.yearFraction(evaluationDate, from);
		//	Real to_t = dc.yearFraction(evaluationDate, to);

		//	Real r1 = this->rate(path, from_t);
		//	Real r2 = this->rate(path, to_t);
		//	Real r = 0.5 * (r1 + r2);
		//	return exp(r*fraction);

		//}

		Real compound(const Path& path, Date from, Date to, DayCounter dc) const
		{
			Date evaluationDate = Settings::instance().evaluationDate();

			Real fraction = dc.yearFraction(from, to);
			Real from_t = dc.yearFraction(evaluationDate, from);
			Real to_t = dc.yearFraction(evaluationDate, to);

			const TimeGrid& timeGrid = path.timeGrid();

			Size from_pos = timeGrid.closestIndex(from_t, -1);
			Size to_pos = timeGrid.closestIndex(to_t, -1);

			if (from_pos == to_pos)
				return this->rate(path, from_t);

			Real sum = 0.0;
			for (Size i = from_pos; i < to_pos; i++)
			{
				Real r = path[i];
				sum += r;
			}

			sum = sum / (to_pos - from_pos);

			return exp(sum * fraction);

		}

		Real discount(const Path& path, Date from, Date to, DayCounter dc) const
		{
			return 1.0/compound(path, from, to, dc);

		}

		// t 기준의 amount // float 부리 해야함.
		//Real legFltAmt(const Leg& leg, const Path& path, Date excerciseDate, Date next_excerciseDate) const
		//{
		//	Date evaluationDate = Settings::instance().evaluationDate();

		//	std::vector<Size> ex_t_pos;
		//	std::vector<Real> amounts;
		//	std::vector<Real> rates;

		//	Size legNum = leg.size();
		//	Real amount = 0.0;
		//	Size startLegCount = 0;

		//	for (Size i = 0; i < legNum; i++)
		//	{
		//		startLegCount = i;
		//		Date d = leg[i]->date();
		//		if (evaluationDate < d) // 당일꺼는 안넣음
		//			break;
		//	}

		//	//for (Size j = startLegCount; j < legNum; j++)
		//	for (Size j = startLegCount; j < legNum; j++)
		//	{
		//		boost::shared_ptr<FloatingRateCoupon> fltCpn
		//			= boost::dynamic_pointer_cast<FloatingRateCoupon>(leg[j]);
		//		Date paymentDate = fltCpn->date();
		//
		//		// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
		//		if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
		//			continue;

		//		DayCounter dc = fltCpn->dayCounter();
		//
		//		Real t = dc.yearFraction(evaluationDate, paymentDate);
		//		Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

		//		if (reset_t < 0)
		//		{
		//			amount += fltCpn->amount();
		//		}
		//		else
		//		{
		//			Real notional = fltCpn->notional();
		//			Size pos = path.timeGrid().closestIndex(reset_t);
		//			Real rate = this->rate(path, reset_t);
		//			Real accrualPeriod = fltCpn->accrualPeriod();
		//			Real compound = 1.0 / this->hw_model_->discountBond(reset_t, reset_t + accrualPeriod, rate);
		//			Real float_rate = std::pow(compound, 1.0 / accrualPeriod)- 1.0;
		//
		//			Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
		//			Real coupon_discounted = coupon * this->discount(path, excerciseDate, paymentDate, dc); // 할인함
		//			amount += coupon_discounted;
		//		}

		//	}

		//	return amount;
		//}

		Real legFltAmt(const Leg& leg, const Path& path, Date excerciseDate, Date next_excerciseDate, Size& couponCalculatedCount) const
		{
			couponCalculatedCount = 0;

			Date evaluationDate = Settings::instance().evaluationDate();

			std::vector<Size> ex_t_pos;
			std::vector<Real> amounts;
			std::vector<Real> rates;

			Size legNum = leg.size();
			Real amount = 0.0;
			Size startLegCount = 0;

			for (Size i = 0; i < legNum; i++)
			{
				startLegCount = i;
				Date d = leg[i]->date();
				if (evaluationDate < d) // 당일꺼는 안넣음
					break;
			}

			//for (Size j = startLegCount; j < legNum; j++)
			for (Size j = startLegCount; j < legNum; j++)
			{
				boost::shared_ptr<FloatingRateCoupon> fltCpn
					= boost::dynamic_pointer_cast<FloatingRateCoupon>(leg[j]);
				Date paymentDate = fltCpn->date();

				// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
				if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
					continue;

				couponCalculatedCount += 1;
				DayCounter dc = fltCpn->dayCounter();

				Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
				Real payment_t = dc.yearFraction(evaluationDate, paymentDate);
				Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

				Real rate = this->rate(path, ex_t);

				if (reset_t < 0)
				{
					//Real accrualPeriod = fltCpn->accrualPeriod() - fltCpn->accruedPeriod(evaluationDate);
					//amount += fltCpn->nominal() * fltCpn->rate() * accrualPeriod * this->hw_model_->discountBond(ex_t, payment_t, rate);
					amount += fltCpn->amount() * this->hw_model_->discountBond(ex_t, payment_t, rate);
				}
				else if (reset_t < ex_t)
				{
					Real notional = fltCpn->notional();
					Size pos = path.timeGrid().closestIndex(reset_t);
					Real reset_rate = this->rate(path, reset_t);
					Real accrualPeriod = fltCpn->accrualPeriod();
					Real compound = 1.0 / this->hw_model_->discountBond(reset_t, reset_t + accrualPeriod, reset_rate);
					//Real float_rate = std::pow(compound, 1.0 / accrualPeriod)- 1.0;
					Real float_rate = (compound - 1.0)/accrualPeriod;

					Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
					Real coupon_discounted = coupon * this->hw_model_->discountBond(ex_t, payment_t, rate); // 할인함
					amount += coupon_discounted;
				}
				else
				{
					Real notional = fltCpn->notional();

					Real accrualPeriod = fltCpn->accrualPeriod();
					Real discount_reset = this->hw_model_->discountBond(ex_t, reset_t, rate);
					Real discount_payment = this->hw_model_->discountBond(ex_t, payment_t, rate);
					Real compound = discount_reset / discount_payment;
					//Real float_rate = std::pow(compound, 1.0 / accrualPeriod) - 1.0;
					Real float_rate = (compound - 1.0) / accrualPeriod;
					Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
					Real coupon_discounted = coupon * discount_payment; // 할인함

					amount += coupon_discounted;
				}

			}

			return amount;
		}

		//Real legFixAmt(const Leg& leg, const Path& path, Date excerciseDate, Size legCount) const
		//{
		//	DayCounter dc = Thirty360();

		//	boost::shared_ptr<SimpleCashFlow> simpleCpn
		//		= boost::dynamic_pointer_cast<SimpleCashFlow>(leg[legCount]);

		//	Real df = this->discount(path, excerciseDate, simpleCpn->date(), dc);

		//	return simpleCpn->amount() * df;
		//}


		// t 기준의 amount
		Real legFixAmt(const Leg& leg, const Path& path, Date excerciseDate, Size legCount) const
		{
			Date evaluationDate = Settings::instance().evaluationDate();

			DayCounter dc = Thirty360();

			boost::shared_ptr<SimpleCashFlow> simpleCpn
				= boost::dynamic_pointer_cast<SimpleCashFlow>(leg[legCount]);

			Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
			Real payment_t = dc.yearFraction(evaluationDate, simpleCpn->date());
			Real rate = this->rate(path, ex_t);

			Real df = this->hw_model_->discountBond(ex_t, payment_t, rate);

			return simpleCpn->amount() * df;
		}

		std::vector<Real> discounts(const Path& path, std::vector<Date> excerciseDates) const
		{
			Date evaluationDate = Settings::instance().evaluationDate();

			std::vector<Real> res;
			Size len = excerciseDates.size();

			TimeGrid timeGrid = path.timeGrid();
			DayCounter dc = Actual365Fixed();

			for (Size i = 0; i < len; i++)
			{
				Date date = excerciseDates[i];
				if (date < evaluationDate)
					continue;

				Real sum = 0.0;
				Real t = dc.yearFraction(evaluationDate, date);
				Size pos = timeGrid.closestIndex(t, -1);
				Real pos_t = timeGrid[pos];

				Real check_accural = 0.0;

				for (Size j = 0; j < pos; j++)
				{
					sum += path[j] * timeGrid.dt(j);
					check_accural += timeGrid.dt(j);
				}

				sum += (path[pos] * (t - timeGrid[pos]));
				check_accural += (t - timeGrid[pos]);

				Real diff = t - check_accural;

				res.push_back(exp(-sum));
			}

			return res;

		}

		void calculate() const
		{
			Real side = 1.0;

			Date evaluationDate = Settings::instance().evaluationDate();

			std::vector<boost::function1<Real, Real>> v =
				LsmBasisSystem::pathBasisSystem(2, LsmBasisSystem::Hermite);

			const Size n = this->samples_;
			std::vector<Real> prices(n), exercise(n);

			// precalculate!
			Leg fixedLeg = this->arguments_.fixedLeg;
			Leg floatingLeg = this->arguments_.floatingLeg;

			DayCounter dc = Actual365Fixed();
			DayCounter dc_360 = Thirty360();
			Time maturityTime = dc.yearFraction(evaluationDate, this->arguments_.excerciseDates.back());

			TimeGrid timeGrid(maturityTime, this->timeSteps_);
			this->process_->fitting()->set_timeGrid(timeGrid);

			PseudoRandom p_rand = PseudoRandom();

			PseudoRandom::rsg_type gen =
				p_rand.make_sequence_generator(timeGrid.size()-1, this->samples_);

			SingleVariate<PseudoRandom>::path_generator_type pathGen
				= SingleVariate<PseudoRandom>::path_generator_type(this->process_ ,
					timeGrid,
					gen,
					false);

			// std::vector<Time> excerciseTimes;
			const std::vector<Date>& excerciseDates = this->arguments_.excerciseDates;
			const Size len = excerciseDates.size();
			coeff_.reset(new Array[len-1]);

			std::vector<std::vector<Real>> precalculated_fixed_amts;
			std::vector<std::vector<Real>> precalculated_float_amts;
			std::vector<std::vector<Real>> precalculated_discounts;

			std::vector<Path> paths;

			std::vector<Real> x;
			std::vector<Real> y;

			// zero callable은 마지막에서 땡기면서 올때, floating cf가 계속 빠지는 구조
			// 처음에 par로거래하고 floating을 주면서 계속 npv가 올라가는 구조여서 그럼.
			for (Size j = 0; j < n; ++j)
			{
				paths.push_back(pathGen.next().value);
				prices[j] = - fixedLeg[len - 1]->amount();

				 // std::cout << j << std::endl;
			}

			std::cout << "--------------------------------------------------" << std::endl;
			Size couponCalculatedCount = 0;

			for (int i = len - 2; i >= 0; --i) {
				y.clear();
				x.clear();

				// std::cout << i << std::endl;

				Time t = dc.yearFraction(evaluationDate, excerciseDates[i]);
				Time next_t = dc.yearFraction(evaluationDate, excerciseDates[i+1]);
				//this->hw_model_ = boost::shared_ptr<HullWhite>(
				//	new HullWhite(this->process_->fitting(), this->process_->a(t), this->process_->sigma(t)));

				//roll back step
				for (Size j = 0; j < n; ++j)
				{
					Path path = paths[j];

					Real rate = this->rate(path, t);
					//Real df = this->discount(path, excerciseDates[i], excerciseDates[i + 1], dc_360);
					Real df = this->hw_model_->discountBond(t, next_t, rate);

					// discounted(excerciseDates[i])
					// Real fixed_amt = this->legFixAmt(fixedLeg, path, excerciseDates[i], i+1);
					Real float_amt = this->legFltAmt(floatingLeg, path, excerciseDates[i], excerciseDates[i+1], couponCalculatedCount);

					Real price = prices[j] * df + float_amt;

					x.push_back(rate);
					y.push_back(price);
				}

				coeff_[i] = GeneralLinearLeastSquares(x, y, v).coefficients();

				for (Size j = 0; j < n; ++j) {
					Path path = paths[j];
					Real rate = this->rate(path, t);

					// Real fixed_amt = this->legFixAmt(fixedLeg, path, excerciseDates[i], i + 1);
					Real float_amt = this->legFltAmt(floatingLeg, path, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);
					//Real df = this->discount(path, excerciseDates[i], excerciseDates[i + 1], dc_360);
					Real df = this->hw_model_->discountBond(t, next_t, rate);

					Real price = prices[j] * df + float_amt;
					prices[j] = price;

					exercise[j] = - fixedLeg[i]->amount();

					Real continuationValue = 0.0;

					for (Size l = 0; l < v.size(); ++l) {
						continuationValue += coeff_[i][l] * v[l](rate);
					}

					if (continuationValue < exercise[j]) {
						prices[j] = exercise[j];
					}
				}
			}

			// 평가일까지 한번 더 땡겨야함.
			Path path = paths[0]; // 어차피 첫번째꺼 사용할 거임. 같은값나오므로.
			Real float_amt = this->legFltAmt(floatingLeg, path, evaluationDate, excerciseDates[0], couponCalculatedCount);
			Real df = this->process_->fitting()->discount(excerciseDates[0]);

			for (Size j = 0; j < n; ++j) {
				// Real df = this->discount(path, evaluationDate, excerciseDates[0], dc);

				Real price = prices[j] * df + float_amt;
				prices[j] = price;
			}

			Real npv = 0.0;

			for (Size i = 0; i < n; i++)
			{
				npv += prices[i] / n;
			}

			this->results_.value = npv ;

			std::cout << npv << std::endl;

			//Time maturity = Thirty360().yearFraction(evaluationDate, fixedLeg.back()->date());
			//Real df = process_->fitting()->discount(maturity);
			//std::cout << fixedLeg.back()->amount() * df << " - " << fixedLeg.back()->amount() << " - " << df << std::endl;

			//std::cout << CashFlows::npv(floatingLeg, *(process_->fitting().currentLink()), false) << std::endl;
		}

	private:
		Size samples_;
		Size polynomOrder_;
		// LsmBasisSystem::PolynomType polynomType_;
		mutable boost::scoped_array<Array> coeff_;

		boost::shared_ptr<HullWhite1FactorProcess> process_;
		mutable boost::shared_ptr<HullWhite> hw_model_;

		Size timeSteps_;
		Size seed_;
		Size nCalibrationSamples_;

	};


	class MCZeroCouponCallableSwapEngine : public TestEngine::ZeroCouponCallableSwap::engine {
	public:
		MCZeroCouponCallableSwapEngine(
			//const boost::shared_ptr<StochasticProcess>& process,
			const boost::shared_ptr<HullWhite1FactorProcess>& process,
			Size timeSteps,
			Size samples,
			BigNatural seed,
			Size nCalibrationSamples)
			: process_(process),
			timeSteps_(timeSteps), samples_(samples), seed_(seed), nCalibrationSamples_(nCalibrationSamples)
		{
			this->hw_model_ = boost::shared_ptr<HullWhite>(
				new HullWhite(process->fitting(), process->a(10.0), process->sigma(10.0)));
		}

		Real legNpv(const Leg& leg, const Path& path, Time t)
		{
			return 0.0;
		}

		Real rate(const Path& path, Time t) const
		{
			Size pos = path.timeGrid().closestIndex(t, -1);

			Real r1 = path[pos];
			Real r2 = path[pos + 1];
			Real t1 = path.timeGrid().at(pos);
			Real t2 = path.timeGrid().at(pos + 1);
			Real res = r1 * (t2 - t) + r2 * (t - t1);

			return res / (t2 - t1);

		}

		Real compound(const Path& path, Date from, Date to, DayCounter dc) const
		{
			Date evaluationDate = Settings::instance().evaluationDate();

			Real fraction = dc.yearFraction(from, to);
			Real from_t = dc.yearFraction(evaluationDate, from);
			Real to_t = dc.yearFraction(evaluationDate, to);

			const TimeGrid& timeGrid = path.timeGrid();

			Size from_pos = timeGrid.closestIndex(from_t, -1);
			Size to_pos = timeGrid.closestIndex(to_t, -1);

			if (from_pos == to_pos)
				return this->rate(path, from_t);

			Real sum = 0.0;
			for (Size i = from_pos; i < to_pos; i++)
			{
				Real r = path[i];
				sum += r;
			}

			sum = sum / (to_pos - from_pos);

			return exp(sum * fraction);

		}

		Real discount(const Path& path, Date from, Date to, DayCounter dc) const
		{
			return 1.0 / compound(path, from, to, dc);

		}

		// t 기준의 amount // float 부리 해야함.
		//Real legFltAmt(const Leg& leg, const Path& path, Date excerciseDate, Date next_excerciseDate) const
		//{
		//	Date evaluationDate = Settings::instance().evaluationDate();

		//	std::vector<Size> ex_t_pos;
		//	std::vector<Real> amounts;
		//	std::vector<Real> rates;

		//	Size legNum = leg.size();
		//	Real amount = 0.0;
		//	Size startLegCount = 0;

		//	for (Size i = 0; i < legNum; i++)
		//	{
		//		startLegCount = i;
		//		Date d = leg[i]->date();
		//		if (evaluationDate < d) // 당일꺼는 안넣음
		//			break;
		//	}

		//	//for (Size j = startLegCount; j < legNum; j++)
		//	for (Size j = startLegCount; j < legNum; j++)
		//	{
		//		boost::shared_ptr<FloatingRateCoupon> fltCpn
		//			= boost::dynamic_pointer_cast<FloatingRateCoupon>(leg[j]);
		//		Date paymentDate = fltCpn->date();
		//
		//		// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
		//		if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
		//			continue;

		//		DayCounter dc = fltCpn->dayCounter();
		//
		//		Real t = dc.yearFraction(evaluationDate, paymentDate);
		//		Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

		//		if (reset_t < 0)
		//		{
		//			amount += fltCpn->amount();
		//		}
		//		else
		//		{
		//			Real notional = fltCpn->notional();
		//			Size pos = path.timeGrid().closestIndex(reset_t);
		//			Real rate = this->rate(path, reset_t);
		//			Real accrualPeriod = fltCpn->accrualPeriod();
		//			Real compound = 1.0 / this->hw_model_->discountBond(reset_t, reset_t + accrualPeriod, rate);
		//			Real float_rate = std::pow(compound, 1.0 / accrualPeriod)- 1.0;
		//
		//			Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
		//			Real coupon_discounted = coupon * this->discount(path, excerciseDate, paymentDate, dc); // 할인함
		//			amount += coupon_discounted;
		//		}

		//	}

		//	return amount;
		//}

		Real legFltAmt(const Leg& leg, const Path& path, Date excerciseDate, Date next_excerciseDate, Size& couponCalculatedCount) const
		{
			couponCalculatedCount = 0;

			Date evaluationDate = Settings::instance().evaluationDate();

			std::vector<Size> ex_t_pos;
			std::vector<Real> amounts;
			std::vector<Real> rates;

			Size legNum = leg.size();
			Real amount = 0.0;
			Size startLegCount = 0;

			for (Size i = 0; i < legNum; i++)
			{
				startLegCount = i;
				Date d = leg[i]->date();
				if (evaluationDate < d) // 당일꺼는 안넣음
					break;
			}

			//for (Size j = startLegCount; j < legNum; j++)
			for (Size j = startLegCount; j < legNum; j++)
			{
				boost::shared_ptr<FloatingRateCoupon> fltCpn
					= boost::dynamic_pointer_cast<FloatingRateCoupon>(leg[j]);
				Date paymentDate = fltCpn->date();

				// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
				if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
					continue;

				couponCalculatedCount += 1;
				DayCounter dc = fltCpn->dayCounter();

				Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
				Real payment_t = dc.yearFraction(evaluationDate, paymentDate);
				Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

				Real rate = this->rate(path, ex_t);
				Real df = this->discount(path, excerciseDate, paymentDate, dc);

				if (reset_t < 0)
				{
					amount += fltCpn->amount() * df;
				}
				else
				{
					Real notional = fltCpn->notional();
					Size pos = path.timeGrid().closestIndex(reset_t);
					Real reset_rate = this->rate(path, reset_t);
					Real accrualPeriod = fltCpn->accrualPeriod();
					Real compound = 1.0 / this->hw_model_->discountBond(reset_t, reset_t + accrualPeriod, reset_rate);
					Real float_rate = (compound - 1.0) / accrualPeriod;

					Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
					Real coupon_discounted = coupon * this->discount(path, excerciseDate, paymentDate, dc); // 할인함
					amount += coupon_discounted;
				}
			}

			return amount;
		}

		void calculate() const
		{
			Real side = 1.0;

			Date evaluationDate = Settings::instance().evaluationDate();

			std::vector<boost::function1<Real, Real>> v =
				LsmBasisSystem::pathBasisSystem(2, LsmBasisSystem::Hermite);

			const Size n = this->samples_;
			std::vector<Real> prices(n), exercise(n);

			// precalculate!
			Leg fixedLeg = this->arguments_.fixedLeg;
			Leg floatingLeg = this->arguments_.floatingLeg;

			DayCounter dc = Actual365Fixed();
			DayCounter dc_360 = Thirty360();
			Time maturityTime = dc.yearFraction(evaluationDate, this->arguments_.excerciseDates.back());

			TimeGrid timeGrid(maturityTime, this->timeSteps_);
			this->process_->fitting()->set_timeGrid(timeGrid);

			PseudoRandom p_rand = PseudoRandom();

			PseudoRandom::rsg_type gen =
				p_rand.make_sequence_generator(timeGrid.size() - 1, this->samples_);

			SingleVariate<PseudoRandom>::path_generator_type pathGen
				= SingleVariate<PseudoRandom>::path_generator_type(this->process_,
					timeGrid,
					gen,
					false);

			// std::vector<Time> excerciseTimes;
			const std::vector<Date>& excerciseDates = this->arguments_.excerciseDates;
			const Size len = excerciseDates.size();
			coeff_.reset(new Array[len - 1]);

			std::vector<std::vector<Real>> precalculated_fixed_amts;
			std::vector<std::vector<Real>> precalculated_float_amts;
			std::vector<std::vector<Real>> precalculated_discounts;

			std::vector<Path> paths;

			std::vector<Real> x;
			std::vector<Real> y;

			// zero callable은 마지막에서 땡기면서 올때, floating cf가 계속 빠지는 구조
			// 처음에 par로거래하고 floating을 주면서 계속 npv가 올라가는 구조여서 그럼.
			for (Size j = 0; j < n; ++j)
			{
				paths.push_back(pathGen.next().value);
				prices[j] = -fixedLeg[len - 1]->amount();
			}

			std::cout << "--------------------------------------------------" << std::endl;
			Size couponCalculatedCount = 0;

			for (int i = len - 2; i >= 0; --i) {
				y.clear();
				x.clear();

				Time t = dc.yearFraction(evaluationDate, excerciseDates[i]);
				//this->hw_model_ = boost::shared_ptr<HullWhite>(
				//	new HullWhite(this->process_->fitting(), this->process_->a(t), this->process_->sigma(t)));

				//roll back step
				for (Size j = 0; j < n; ++j)
				{
					Path path = paths[j];

					Real rate = this->rate(path, t);
					Real df = this->discount(path, excerciseDates[i], excerciseDates[i + 1], dc_360);
					Real float_amt = this->legFltAmt(floatingLeg, path, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);

					Real price = prices[j] * df + float_amt;

					x.push_back(rate);
					y.push_back(price);
				}

				coeff_[i] = GeneralLinearLeastSquares(x, y, v).coefficients();

				for (Size j = 0; j < n; ++j) {
					Path path = paths[j];
					Real rate = this->rate(path, t);

					Real float_amt = this->legFltAmt(floatingLeg, path, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);
					Real df = this->discount(path, excerciseDates[i], excerciseDates[i + 1], dc_360);

					Real price = prices[j] * df + float_amt;
					prices[j] = price;

					exercise[j] = -fixedLeg[i]->amount();

					Real continuationValue = 0.0;

					for (Size l = 0; l < v.size(); ++l) {
						continuationValue += coeff_[i][l] * v[l](rate);
					}

					if (continuationValue < exercise[j]) {
						prices[j] = exercise[j];
					}
				}
			}

			// 평가일까지 한번 더 땡겨야함.
			Path path = paths[0]; // 어차피 첫번째꺼 사용할 거임. 같은값나오므로.
			Real float_amt = this->legFltAmt(floatingLeg, path, evaluationDate, excerciseDates[0], couponCalculatedCount);
			Real df = this->discount(path, evaluationDate, excerciseDates[0], dc);

			for (Size j = 0; j < n; ++j) {

				Real price = prices[j] * df + float_amt;
				prices[j] = price;
			}

			Real npv = 0.0;

			for (Size i = 0; i < n; i++)
			{
				npv += prices[i] / n;
			}

			this->results_.value = npv;

			std::cout << npv << std::endl;

			//Time maturity = Thirty360().yearFraction(evaluationDate, fixedLeg.back()->date());
			//Real df = process_->fitting()->discount(maturity);
			//std::cout << fixedLeg.back()->amount() * df << " - " << fixedLeg.back()->amount() << " - " << df << std::endl;

			//std::cout << CashFlows::npv(floatingLeg, *(process_->fitting().currentLink()), false) << std::endl;
		}

	private:
		Size samples_;
		Size polynomOrder_;
		// LsmBasisSystem::PolynomType polynomType_;
		mutable boost::scoped_array<Array> coeff_;

		boost::shared_ptr<HullWhite1FactorProcess> process_;
		mutable boost::shared_ptr<HullWhite> hw_model_;

		Size timeSteps_;
		Size seed_;
		Size nCalibrationSamples_;

	};

}

namespace TestEngine
{
	Size simulNum()
	{
		return 10000;
	}

	boost::shared_ptr<YieldTermStructure> irs_curve()
	{
		Date evaluationDate = Settings::instance().evaluationDate();
		std::vector<MarketCurveRate> marketCurveRates
		{
			MarketCurveRate("1D", 0.0049, MarketCurveRate::Type::Cash),
			MarketCurveRate("3M", 0.0073, MarketCurveRate::Type::Cash),
			MarketCurveRate("6M", 0.007225, MarketCurveRate::Type::Swap),
			MarketCurveRate("9M", 0.007575, MarketCurveRate::Type::Swap),
			MarketCurveRate("1Y", 0.008125, MarketCurveRate::Type::Swap),
			MarketCurveRate("18M", 0.009425, MarketCurveRate::Type::Swap),
			MarketCurveRate("2Y", 0.01045, MarketCurveRate::Type::Swap),
			MarketCurveRate("3Y", 0.012425, MarketCurveRate::Type::Swap),
			MarketCurveRate("4Y", 0.014, MarketCurveRate::Type::Swap),
			MarketCurveRate("5Y", 0.015, MarketCurveRate::Type::Swap),
			MarketCurveRate("6Y", 0.01565, MarketCurveRate::Type::Swap),
			MarketCurveRate("7Y", 0.016175, MarketCurveRate::Type::Swap),
			MarketCurveRate("8Y", 0.016525, MarketCurveRate::Type::Swap),
			MarketCurveRate("9Y", 0.016975, MarketCurveRate::Type::Swap),
			MarketCurveRate("10Y", 0.017275, MarketCurveRate::Type::Swap),
			MarketCurveRate("12Y", 0.017725, MarketCurveRate::Type::Swap),
			MarketCurveRate("15Y", 0.0176, MarketCurveRate::Type::Swap),
			MarketCurveRate("20Y", 0.017525, MarketCurveRate::Type::Swap),
			MarketCurveRate("25Y", 0.01765, MarketCurveRate::Type::Swap),
			MarketCurveRate("30Y", 0.01765, MarketCurveRate::Type::Swap)
		};

		return YieldCurveExt::bootstrapping_ccp(evaluationDate, marketCurveRates,
			Interpolator1D::Linear,
			Extrapolator1D::FlatForward, "IRSKRW_KRCCP");

	}

	boost::shared_ptr<HullWhiteProcess> process_hw_const()
	{
		Real alpha = 0.05;
		Real sigma = 0.006188;
		boost::shared_ptr<HullWhiteProcess> process(
			new HullWhiteProcess(Handle<YieldTermStructure>(irs_curve()), alpha, sigma));

		return process;

	}

	//boost::shared_ptr<HullWhite1FactorProcess> process_positive_hw_timevaring()
	boost::shared_ptr<HullWhite1FactorProcess> process()
	{
		std::vector<Real> alpha_times{ 100.0 };
		PiecewiseConstantParameter alpha_para(alpha_times);
		alpha_para.setParam(0, 0.06188);

		std::vector<Real> sigma_times{ 0.083333,0.25,0.5,0.75,1,1.5,2,3,4,5,7,10 };

		// 이거는 simple하게 구한 볼
		std::vector<Real> sigma_values
		{
			0.006188,
			0.00580129905905104,
			0.00566797309450212,
			0.00539632208082505,
			0.00509562518244817,
			0.00490434858059661,
			0.00438993086505927,
			0.00461860790715125,
			0.00428171893052311,
			0.00393527394217989,
			0.004635,
			0.004635
		};

		// 이거는 복잡하게 구한 포워드 볼
		//std::vector<Real> sigma_values
		//{
		//	0.006188,
		//	0.00577550985818638,
		//	0.00559417299802497,
		//	0.00525259554356766,
		//	0.00488545691359162,
		//	0.00463191590425029,
		//	0.00399507392210221,
		//	0.00412757749289384,
		//	0.00357839278576286,
		//	0.00303426929987496,
		//	0.00360974162953596,
		//	0.00326622929584623
		//};

		// 이건 spot vol
		//std::vector<Real> sigma_values
		//{
		//	0.006188,
		//	0.005933,
		//	0.005802,
		//	0.00567,
		//	0.005532,
		//	0.005331,
		//	0.005112,
		//	0.004953,
		//	0.004794,
		//	0.004635,
		//	0.004635,
		//	0.004635
		//};

		LinearInterpolationParameter sigma_para(sigma_times, sigma_values);

		//boost::shared_ptr<HullWhite1FactorProcess> process(
		//	new HullWhite1FactorProcess(Handle<YieldTermStructure>(irs_curve()), alpha_para, sigma_para));
		boost::shared_ptr<HullWhite1FactorProcess> process(
			new HullWhite1FactorPositiveProcess(Handle<YieldTermStructure>(irs_curve()), alpha_para, sigma_para));

		return process;

	}

	boost::shared_ptr<HullWhite1FactorProcess> process_const2()
	{
		std::vector<Real> times{100.0};
		PiecewiseConstantParameter alpha_para(times);
		alpha_para.setParam(0, 0.05);

		PiecewiseConstantParameter sigma_para(times);
		sigma_para.setParam(0, 0.006188);

		Handle<YieldTermStructure> fitting_h(irs_curve());

		boost::shared_ptr<HullWhite1FactorProcess> process(
			new HullWhite1FactorProcess(fitting_h, alpha_para, sigma_para));

		return process;

	}

	//boost::shared_ptr<HullWhite1FactorProcess> process()
	boost::shared_ptr<HullWhite1FactorProcess> process_const_positive()
	{
		std::vector<Real> times{ 100.0 };
		PiecewiseConstantParameter alpha_para(times);
		alpha_para.setParam(0, 0.05);

		PiecewiseConstantParameter sigma_para(times);
		sigma_para.setParam(0, 0.006188);

		Handle<YieldTermStructure> fitting_h(irs_curve());

		boost::shared_ptr<HullWhite1FactorProcess> process(
			new HullWhite1FactorPositiveProcess(fitting_h, alpha_para, sigma_para));

		return process;

	}

	void test_28888()
	{
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-04-30");
		Date evaluationDate = Settings::instance().evaluationDate();
		// get curve 랑 등등 해서 완성했어.
		// 이거는 여기다가 그냥 hard coding 해서 넣는 방식임. 그래서 여러가지 바꿔보는 것...?
		Real notional = 60000000000;

		Leg fixedLeg;
		Real rate = 0.0212;
		std::vector<Date> fixed_payment_dates
		{
			//DateParser::parseISO("2017-08-23"),
			//DateParser::parseISO("2018-08-23"),
			//DateParser::parseISO("2019-08-23"),
			//DateParser::parseISO("2020-08-24"),
			DateParser::parseISO("2021-07-20"),
			DateParser::parseISO("2022-07-20"),
			DateParser::parseISO("2023-07-20"),
			DateParser::parseISO("2024-07-22"),
			DateParser::parseISO("2025-07-21"),
			DateParser::parseISO("2026-07-20"),
			DateParser::parseISO("2027-07-20"),
			DateParser::parseISO("2028-07-20"),
			DateParser::parseISO("2029-07-20"),
			DateParser::parseISO("2030-07-22"),
			DateParser::parseISO("2031-07-21"),
			DateParser::parseISO("2032-07-20"),
			DateParser::parseISO("2033-07-20"),
			DateParser::parseISO("2034-07-20"),
			DateParser::parseISO("2035-07-20"),
			DateParser::parseISO("2036-07-21")
		};

		for (Size i = 0; i < fixed_payment_dates.size(); i++)
		{
			Real amount = notional * rate * (i + 5);
			boost::shared_ptr<SimpleCashFlow> scf(new SimpleCashFlow(amount, fixed_payment_dates[i]));
			fixedLeg.push_back(scf);
		}

		Leg floatingLeg;

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new Libor("KRWCD3M", Period(3, TimeUnit::Months), 1,
				KRWCurrency(), SouthKorea(), Actual365Fixed(), yts_h));

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = -0.0006;

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2016-10-20"), notional, DateParser::parseISO("2016-07-20"), DateParser::parseISO("2016-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-01-20"), notional, DateParser::parseISO("2016-10-20"), DateParser::parseISO("2017-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-04-20"), notional, DateParser::parseISO("2017-01-20"), DateParser::parseISO("2017-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-07-20"), notional, DateParser::parseISO("2017-04-20"), DateParser::parseISO("2017-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-10-20"), notional, DateParser::parseISO("2017-07-20"), DateParser::parseISO("2017-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-01-22"), notional, DateParser::parseISO("2017-10-20"), DateParser::parseISO("2018-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-04-20"), notional, DateParser::parseISO("2018-01-22"), DateParser::parseISO("2018-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-07-20"), notional, DateParser::parseISO("2018-04-20"), DateParser::parseISO("2018-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-10-22"), notional, DateParser::parseISO("2018-07-20"), DateParser::parseISO("2018-10-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-01-21"), notional, DateParser::parseISO("2018-10-22"), DateParser::parseISO("2019-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-04-22"), notional, DateParser::parseISO("2019-01-21"), DateParser::parseISO("2019-04-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-07-22"), notional, DateParser::parseISO("2019-04-22"), DateParser::parseISO("2019-07-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-10-21"), notional, DateParser::parseISO("2019-07-22"), DateParser::parseISO("2019-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-01-20"), notional, DateParser::parseISO("2019-10-21"), DateParser::parseISO("2020-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-04-20"), notional, DateParser::parseISO("2020-01-20"), DateParser::parseISO("2020-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-07-20"), notional, DateParser::parseISO("2020-04-20"), DateParser::parseISO("2020-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-10-20"), notional, DateParser::parseISO("2020-07-20"), DateParser::parseISO("2020-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-01-20"), notional, DateParser::parseISO("2020-10-20"), DateParser::parseISO("2021-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-04-20"), notional, DateParser::parseISO("2021-01-20"), DateParser::parseISO("2021-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-07-20"), notional, DateParser::parseISO("2021-04-20"), DateParser::parseISO("2021-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-10-20"), notional, DateParser::parseISO("2021-07-20"), DateParser::parseISO("2021-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-01-20"), notional, DateParser::parseISO("2021-10-20"), DateParser::parseISO("2022-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-04-20"), notional, DateParser::parseISO("2022-01-20"), DateParser::parseISO("2022-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-07-20"), notional, DateParser::parseISO("2022-04-20"), DateParser::parseISO("2022-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-10-20"), notional, DateParser::parseISO("2022-07-20"), DateParser::parseISO("2022-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-01-20"), notional, DateParser::parseISO("2022-10-20"), DateParser::parseISO("2023-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-04-20"), notional, DateParser::parseISO("2023-01-20"), DateParser::parseISO("2023-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-07-20"), notional, DateParser::parseISO("2023-04-20"), DateParser::parseISO("2023-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-10-20"), notional, DateParser::parseISO("2023-07-20"), DateParser::parseISO("2023-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-01-22"), notional, DateParser::parseISO("2023-10-20"), DateParser::parseISO("2024-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-04-22"), notional, DateParser::parseISO("2024-01-22"), DateParser::parseISO("2024-04-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-07-22"), notional, DateParser::parseISO("2024-04-22"), DateParser::parseISO("2024-07-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-10-21"), notional, DateParser::parseISO("2024-07-22"), DateParser::parseISO("2024-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-01-20"), notional, DateParser::parseISO("2024-10-21"), DateParser::parseISO("2025-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-04-21"), notional, DateParser::parseISO("2025-01-20"), DateParser::parseISO("2025-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-07-21"), notional, DateParser::parseISO("2025-04-21"), DateParser::parseISO("2025-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-10-20"), notional, DateParser::parseISO("2025-07-21"), DateParser::parseISO("2025-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-01-20"), notional, DateParser::parseISO("2025-10-20"), DateParser::parseISO("2026-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-04-20"), notional, DateParser::parseISO("2026-01-20"), DateParser::parseISO("2026-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-07-20"), notional, DateParser::parseISO("2026-04-20"), DateParser::parseISO("2026-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-10-20"), notional, DateParser::parseISO("2026-07-20"), DateParser::parseISO("2026-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-01-20"), notional, DateParser::parseISO("2026-10-20"), DateParser::parseISO("2027-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-04-20"), notional, DateParser::parseISO("2027-01-20"), DateParser::parseISO("2027-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-07-20"), notional, DateParser::parseISO("2027-04-20"), DateParser::parseISO("2027-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-10-20"), notional, DateParser::parseISO("2027-07-20"), DateParser::parseISO("2027-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-01-20"), notional, DateParser::parseISO("2027-10-20"), DateParser::parseISO("2028-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-04-20"), notional, DateParser::parseISO("2028-01-20"), DateParser::parseISO("2028-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-07-20"), notional, DateParser::parseISO("2028-04-20"), DateParser::parseISO("2028-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-10-20"), notional, DateParser::parseISO("2028-07-20"), DateParser::parseISO("2028-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-01-22"), notional, DateParser::parseISO("2028-10-20"), DateParser::parseISO("2029-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-04-20"), notional, DateParser::parseISO("2029-01-22"), DateParser::parseISO("2029-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-07-20"), notional, DateParser::parseISO("2029-04-20"), DateParser::parseISO("2029-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-10-22"), notional, DateParser::parseISO("2029-07-20"), DateParser::parseISO("2029-10-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-01-21"), notional, DateParser::parseISO("2029-10-22"), DateParser::parseISO("2030-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-04-22"), notional, DateParser::parseISO("2030-01-21"), DateParser::parseISO("2030-04-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-07-22"), notional, DateParser::parseISO("2030-04-22"), DateParser::parseISO("2030-07-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-10-21"), notional, DateParser::parseISO("2030-07-22"), DateParser::parseISO("2030-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-01-20"), notional, DateParser::parseISO("2030-10-21"), DateParser::parseISO("2031-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-04-21"), notional, DateParser::parseISO("2031-01-20"), DateParser::parseISO("2031-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-07-21"), notional, DateParser::parseISO("2031-04-21"), DateParser::parseISO("2031-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-10-20"), notional, DateParser::parseISO("2031-07-21"), DateParser::parseISO("2031-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-01-20"), notional, DateParser::parseISO("2031-10-20"), DateParser::parseISO("2032-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-04-20"), notional, DateParser::parseISO("2032-01-20"), DateParser::parseISO("2032-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-07-20"), notional, DateParser::parseISO("2032-04-20"), DateParser::parseISO("2032-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-10-20"), notional, DateParser::parseISO("2032-07-20"), DateParser::parseISO("2032-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-01-20"), notional, DateParser::parseISO("2032-10-20"), DateParser::parseISO("2033-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-04-20"), notional, DateParser::parseISO("2033-01-20"), DateParser::parseISO("2033-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-07-20"), notional, DateParser::parseISO("2033-04-20"), DateParser::parseISO("2033-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-10-20"), notional, DateParser::parseISO("2033-07-20"), DateParser::parseISO("2033-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-01-20"), notional, DateParser::parseISO("2033-10-20"), DateParser::parseISO("2034-01-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-04-20"), notional, DateParser::parseISO("2034-01-20"), DateParser::parseISO("2034-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-07-20"), notional, DateParser::parseISO("2034-04-20"), DateParser::parseISO("2034-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-10-20"), notional, DateParser::parseISO("2034-07-20"), DateParser::parseISO("2034-10-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-01-22"), notional, DateParser::parseISO("2034-10-20"), DateParser::parseISO("2035-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-04-20"), notional, DateParser::parseISO("2035-01-22"), DateParser::parseISO("2035-04-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-07-20"), notional, DateParser::parseISO("2035-04-20"), DateParser::parseISO("2035-07-20"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-10-22"), notional, DateParser::parseISO("2035-07-20"), DateParser::parseISO("2035-10-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-01-21"), notional, DateParser::parseISO("2035-10-22"), DateParser::parseISO("2036-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-04-21"), notional, DateParser::parseISO("2036-01-21"), DateParser::parseISO("2036-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-07-21"), notional, DateParser::parseISO("2036-04-21"), DateParser::parseISO("2036-07-21"), fixingDays, index, gearing, spread)));


		std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2021-07-20"),
			DateParser::parseISO("2022-07-20"),
			DateParser::parseISO("2023-07-20"),
			DateParser::parseISO("2024-07-22"),
			DateParser::parseISO("2025-07-21"),
			DateParser::parseISO("2026-07-20"),
			DateParser::parseISO("2027-07-20"),
			DateParser::parseISO("2028-07-20"),
			DateParser::parseISO("2029-07-20"),
			DateParser::parseISO("2030-07-22"),
			DateParser::parseISO("2031-07-21"),
			DateParser::parseISO("2032-07-20"),
			DateParser::parseISO("2033-07-20"),
			DateParser::parseISO("2034-07-20"),
			DateParser::parseISO("2035-07-20"),
			DateParser::parseISO("2036-07-21")
		};

		index->addFixing(DateParser::parseISO("2021-04-19"), 0.0073);

		boost::shared_ptr<LegExerciseOption> option(new LegExerciseOption(optioin_exdates, fixed_payment_dates));

		// 아 그거는 어떻게 하지...? payoff. is zerocallable payoff
		bool isFixedRec = true;

		boost::shared_ptr<Swap> swap(
			new ZeroCouponCallableSwap(fixedLeg, floatingLeg, option, isFixedRec));

		Size timeSteps = Actual365Fixed().yearFraction(evaluationDate, optioin_exdates.back()) * 365;
		Size samples = simulNum();

		boost::shared_ptr<MCZeroCouponCallableSwapEngine> engine(
			new MCZeroCouponCallableSwapEngine(
				process(),
				timeSteps,
				samples,
				1,
				2048));

		swap->setPricingEngine(engine);

		try
		{
			swap->NPV();
			// a360    : -2,470,225,814
			// calypso : -2,563,845,632
			//


		}
		catch (const std::exception& e)
		{
			std::cout << e.what();
		}
	}

	void test_37958()
	{
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-04-30");
		Date evaluationDate = Settings::instance().evaluationDate();
		// get curve 랑 등등 해서 완성했어.
		// 이거는 여기다가 그냥 hard coding 해서 넣는 방식임. 그래서 여러가지 바꿔보는 것...?
		Real notional = 50000000000;

		Leg fixedLeg;
		Real rate = 0.0215;
		std::vector<Date> fixed_payment_dates
		{
			//DateParser::parseISO("2017-08-23"),
			//DateParser::parseISO("2018-08-23"),
			//DateParser::parseISO("2019-08-23"),
			//DateParser::parseISO("2020-08-24"),
			DateParser::parseISO("2021-08-23"),
			DateParser::parseISO("2022-08-23"),
			DateParser::parseISO("2023-08-23"),
			DateParser::parseISO("2024-08-23"),
			DateParser::parseISO("2025-08-25"),
			DateParser::parseISO("2026-08-24"),
			DateParser::parseISO("2027-08-23"),
			DateParser::parseISO("2028-08-23"),
			DateParser::parseISO("2029-08-23"),
			DateParser::parseISO("2030-08-23"),
			DateParser::parseISO("2031-08-25"),
			DateParser::parseISO("2032-08-23"),
			DateParser::parseISO("2033-08-23"),
			DateParser::parseISO("2034-08-23"),
			DateParser::parseISO("2035-08-23"),
			DateParser::parseISO("2036-08-25")
		};

		for (Size i = 0; i < fixed_payment_dates.size(); i++)
		{
			Real amount = notional * rate * (i+5);
			boost::shared_ptr<SimpleCashFlow> scf(new SimpleCashFlow(amount, fixed_payment_dates[i]));
			fixedLeg.push_back(scf);
		}

		Leg floatingLeg;

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new Libor("KRWCD3M", Period(3,TimeUnit::Months), 1,
					   KRWCurrency(), SouthKorea(), Actual365Fixed(), yts_h));

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = -0.0005;

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2016-11-23"), notional, DateParser::parseISO("2016-08-23"), DateParser::parseISO("2016-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-02-23"), notional, DateParser::parseISO("2016-11-23"), DateParser::parseISO("2017-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-05-23"), notional, DateParser::parseISO("2017-02-23"), DateParser::parseISO("2017-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-08-23"), notional, DateParser::parseISO("2017-05-23"), DateParser::parseISO("2017-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-11-23"), notional, DateParser::parseISO("2017-08-23"), DateParser::parseISO("2017-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-02-23"), notional, DateParser::parseISO("2017-11-23"), DateParser::parseISO("2018-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-05-23"), notional, DateParser::parseISO("2018-02-23"), DateParser::parseISO("2018-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-08-23"), notional, DateParser::parseISO("2018-05-23"), DateParser::parseISO("2018-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-11-23"), notional, DateParser::parseISO("2018-08-23"), DateParser::parseISO("2018-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-02-25"), notional, DateParser::parseISO("2018-11-23"), DateParser::parseISO("2019-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-05-23"), notional, DateParser::parseISO("2019-02-25"), DateParser::parseISO("2019-05-27"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-08-23"), notional, DateParser::parseISO("2019-05-23"), DateParser::parseISO("2019-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-11-25"), notional, DateParser::parseISO("2019-08-23"), DateParser::parseISO("2019-11-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-02-24"), notional, DateParser::parseISO("2019-11-25"), DateParser::parseISO("2020-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-05-25"), notional, DateParser::parseISO("2020-02-24"), DateParser::parseISO("2020-05-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-08-24"), notional, DateParser::parseISO("2020-05-25"), DateParser::parseISO("2020-08-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-11-23"), notional, DateParser::parseISO("2020-08-24"), DateParser::parseISO("2020-11-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-02-23"), notional, DateParser::parseISO("2020-11-23"), DateParser::parseISO("2021-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-05-24"), notional, DateParser::parseISO("2021-02-23"), DateParser::parseISO("2021-05-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-08-23"), notional, DateParser::parseISO("2021-05-24"), DateParser::parseISO("2021-08-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-11-23"), notional, DateParser::parseISO("2021-08-23"), DateParser::parseISO("2021-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-02-23"), notional, DateParser::parseISO("2021-11-23"), DateParser::parseISO("2022-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-05-23"), notional, DateParser::parseISO("2022-02-23"), DateParser::parseISO("2022-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-08-23"), notional, DateParser::parseISO("2022-05-23"), DateParser::parseISO("2022-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-11-23"), notional, DateParser::parseISO("2022-08-23"), DateParser::parseISO("2022-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-02-23"), notional, DateParser::parseISO("2022-11-23"), DateParser::parseISO("2023-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-05-23"), notional, DateParser::parseISO("2023-02-23"), DateParser::parseISO("2023-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-08-23"), notional, DateParser::parseISO("2023-05-23"), DateParser::parseISO("2023-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-11-23"), notional, DateParser::parseISO("2023-08-23"), DateParser::parseISO("2023-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-02-23"), notional, DateParser::parseISO("2023-11-23"), DateParser::parseISO("2024-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-05-23"), notional, DateParser::parseISO("2024-02-23"), DateParser::parseISO("2024-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-08-23"), notional, DateParser::parseISO("2024-05-23"), DateParser::parseISO("2024-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-11-25"), notional, DateParser::parseISO("2024-08-23"), DateParser::parseISO("2024-11-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-02-24"), notional, DateParser::parseISO("2024-11-25"), DateParser::parseISO("2025-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-05-23"), notional, DateParser::parseISO("2025-02-24"), DateParser::parseISO("2025-05-26"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-08-25"), notional, DateParser::parseISO("2025-05-23"), DateParser::parseISO("2025-08-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-11-24"), notional, DateParser::parseISO("2025-08-25"), DateParser::parseISO("2025-11-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-02-23"), notional, DateParser::parseISO("2025-11-24"), DateParser::parseISO("2026-02-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-05-25"), notional, DateParser::parseISO("2026-02-23"), DateParser::parseISO("2026-05-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-08-24"), notional, DateParser::parseISO("2026-05-25"), DateParser::parseISO("2026-08-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-11-23"), notional, DateParser::parseISO("2026-08-24"), DateParser::parseISO("2026-11-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-02-23"), notional, DateParser::parseISO("2026-11-23"), DateParser::parseISO("2027-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-05-24"), notional, DateParser::parseISO("2027-02-23"), DateParser::parseISO("2027-05-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-08-23"), notional, DateParser::parseISO("2027-05-24"), DateParser::parseISO("2027-08-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-11-23"), notional, DateParser::parseISO("2027-08-23"), DateParser::parseISO("2027-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-02-23"), notional, DateParser::parseISO("2027-11-23"), DateParser::parseISO("2028-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-05-23"), notional, DateParser::parseISO("2028-02-23"), DateParser::parseISO("2028-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-08-23"), notional, DateParser::parseISO("2028-05-23"), DateParser::parseISO("2028-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-11-23"), notional, DateParser::parseISO("2028-08-23"), DateParser::parseISO("2028-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-02-23"), notional, DateParser::parseISO("2028-11-23"), DateParser::parseISO("2029-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-05-23"), notional, DateParser::parseISO("2029-02-23"), DateParser::parseISO("2029-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-08-23"), notional, DateParser::parseISO("2029-05-23"), DateParser::parseISO("2029-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-11-23"), notional, DateParser::parseISO("2029-08-23"), DateParser::parseISO("2029-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-02-25"), notional, DateParser::parseISO("2029-11-23"), DateParser::parseISO("2030-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-05-23"), notional, DateParser::parseISO("2030-02-25"), DateParser::parseISO("2030-05-27"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-08-23"), notional, DateParser::parseISO("2030-05-23"), DateParser::parseISO("2030-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-11-25"), notional, DateParser::parseISO("2030-08-23"), DateParser::parseISO("2030-11-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-02-24"), notional, DateParser::parseISO("2030-11-25"), DateParser::parseISO("2031-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-05-23"), notional, DateParser::parseISO("2031-02-24"), DateParser::parseISO("2031-05-26"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-08-25"), notional, DateParser::parseISO("2031-05-23"), DateParser::parseISO("2031-08-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-11-24"), notional, DateParser::parseISO("2031-08-25"), DateParser::parseISO("2031-11-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-02-23"), notional, DateParser::parseISO("2031-11-24"), DateParser::parseISO("2032-02-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-05-24"), notional, DateParser::parseISO("2032-02-23"), DateParser::parseISO("2032-05-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-08-23"), notional, DateParser::parseISO("2032-05-24"), DateParser::parseISO("2032-08-24"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-11-23"), notional, DateParser::parseISO("2032-08-23"), DateParser::parseISO("2032-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-02-23"), notional, DateParser::parseISO("2032-11-23"), DateParser::parseISO("2033-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-05-23"), notional, DateParser::parseISO("2033-02-23"), DateParser::parseISO("2033-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-08-23"), notional, DateParser::parseISO("2033-05-23"), DateParser::parseISO("2033-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-11-23"), notional, DateParser::parseISO("2033-08-23"), DateParser::parseISO("2033-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-02-23"), notional, DateParser::parseISO("2033-11-23"), DateParser::parseISO("2034-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-05-23"), notional, DateParser::parseISO("2034-02-23"), DateParser::parseISO("2034-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-08-23"), notional, DateParser::parseISO("2034-05-23"), DateParser::parseISO("2034-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-11-23"), notional, DateParser::parseISO("2034-08-23"), DateParser::parseISO("2034-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-02-23"), notional, DateParser::parseISO("2034-11-23"), DateParser::parseISO("2035-02-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-05-23"), notional, DateParser::parseISO("2035-02-23"), DateParser::parseISO("2035-05-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-08-23"), notional, DateParser::parseISO("2035-05-23"), DateParser::parseISO("2035-08-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-11-23"), notional, DateParser::parseISO("2035-08-23"), DateParser::parseISO("2035-11-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-02-25"), notional, DateParser::parseISO("2035-11-23"), DateParser::parseISO("2036-02-25"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-05-23"), notional, DateParser::parseISO("2036-02-25"), DateParser::parseISO("2036-05-26"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-08-25"), notional, DateParser::parseISO("2036-05-23"), DateParser::parseISO("2036-08-25"), fixingDays, index, gearing, spread)));

		std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2021-08-23"),
			DateParser::parseISO("2022-08-23"),
			DateParser::parseISO("2023-08-23"),
			DateParser::parseISO("2024-08-23"),
			DateParser::parseISO("2025-08-25"),
			DateParser::parseISO("2026-08-24"),
			DateParser::parseISO("2027-08-23"),
			DateParser::parseISO("2028-08-23"),
			DateParser::parseISO("2029-08-23"),
			DateParser::parseISO("2030-08-23"),
			DateParser::parseISO("2031-08-25"),
			DateParser::parseISO("2032-08-23"),
			DateParser::parseISO("2033-08-23"),
			DateParser::parseISO("2034-08-23"),
			DateParser::parseISO("2035-08-23"),
			DateParser::parseISO("2036-08-25")
		};

		index->addFixing(DateParser::parseISO("2021-02-22"), 0.0074);

		boost::shared_ptr<LegExerciseOption> option(new LegExerciseOption(optioin_exdates, fixed_payment_dates));

		// 아 그거는 어떻게 하지...? payoff. is zerocallable payoff
		bool isFixedRec = true;

		boost::shared_ptr<Swap> swap(
			new ZeroCouponCallableSwap(fixedLeg, floatingLeg, option, isFixedRec));

		Size timeSteps = Actual365Fixed().yearFraction(evaluationDate, optioin_exdates.back()) * 365;
		Size samples = simulNum();

		boost::shared_ptr<MCZeroCouponCallableSwapEngine> engine(
			new MCZeroCouponCallableSwapEngine(
			process(),
			timeSteps,
			samples,
			1,
			2048));

		swap->setPricingEngine(engine);

		try
		{
			swap->NPV();



		}
		catch (const std::exception& e)
		{
			std::cout << e.what();
		}
	}

	void test_53600()
	{
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-04-30");
		Date evaluationDate = Settings::instance().evaluationDate();
		// get curve 랑 등등 해서 완성했어.
		// 이거는 여기다가 그냥 hard coding 해서 넣는 방식임. 그래서 여러가지 바꿔보는 것...?
		Real notional = 50000000000;

		Leg fixedLeg;
		Real rate = 0.0217;
		std::vector<Date> fixed_payment_dates
		{
			//DateParser::parseISO("2017-08-23"),
			//DateParser::parseISO("2018-08-23"),
			//DateParser::parseISO("2019-08-23"),
			//DateParser::parseISO("2020-08-24"),
			DateParser::parseISO("2021-09-06"),
			DateParser::parseISO("2022-09-05"),
			DateParser::parseISO("2023-09-05"),
			DateParser::parseISO("2024-09-05"),
			DateParser::parseISO("2025-09-05"),
			DateParser::parseISO("2026-09-07"),
			DateParser::parseISO("2027-09-06"),
			DateParser::parseISO("2028-09-05"),
			DateParser::parseISO("2029-09-05"),
			DateParser::parseISO("2030-09-05"),
			DateParser::parseISO("2031-09-05"),
			DateParser::parseISO("2032-09-06"),
			DateParser::parseISO("2033-09-05"),
			DateParser::parseISO("2034-09-05"),
			DateParser::parseISO("2035-09-05"),
			DateParser::parseISO("2036-09-05")
		};

		for (Size i = 0; i < fixed_payment_dates.size(); i++)
		{
			Real amount = notional * rate * (i + 5);
			boost::shared_ptr<SimpleCashFlow> scf(new SimpleCashFlow(amount, fixed_payment_dates[i]));
			fixedLeg.push_back(scf);
		}

		Leg floatingLeg;

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new Libor("KRWCD3M", Period(3, TimeUnit::Months), 1,
				KRWCurrency(), SouthKorea(), Actual365Fixed(), yts_h));

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = -0.0006;

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2016-12-05"), notional, DateParser::parseISO("2016-09-05"), DateParser::parseISO("2016-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-03-06"), notional, DateParser::parseISO("2016-12-05"), DateParser::parseISO("2017-03-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-06-05"), notional, DateParser::parseISO("2017-03-06"), DateParser::parseISO("2017-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-09-05"), notional, DateParser::parseISO("2017-06-05"), DateParser::parseISO("2017-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-12-05"), notional, DateParser::parseISO("2017-09-05"), DateParser::parseISO("2017-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-03-05"), notional, DateParser::parseISO("2017-12-05"), DateParser::parseISO("2018-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-06-05"), notional, DateParser::parseISO("2018-03-05"), DateParser::parseISO("2018-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-09-05"), notional, DateParser::parseISO("2018-06-05"), DateParser::parseISO("2018-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-12-05"), notional, DateParser::parseISO("2018-09-05"), DateParser::parseISO("2018-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-03-05"), notional, DateParser::parseISO("2018-12-05"), DateParser::parseISO("2019-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-06-05"), notional, DateParser::parseISO("2019-03-05"), DateParser::parseISO("2019-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-09-05"), notional, DateParser::parseISO("2019-06-05"), DateParser::parseISO("2019-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-12-05"), notional, DateParser::parseISO("2019-09-05"), DateParser::parseISO("2019-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-03-05"), notional, DateParser::parseISO("2019-12-05"), DateParser::parseISO("2020-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-06-05"), notional, DateParser::parseISO("2020-03-05"), DateParser::parseISO("2020-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-09-07"), notional, DateParser::parseISO("2020-06-05"), DateParser::parseISO("2020-09-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-12-07"), notional, DateParser::parseISO("2020-09-07"), DateParser::parseISO("2020-12-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-03-05"), notional, DateParser::parseISO("2020-12-07"), DateParser::parseISO("2021-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-06-07"), notional, DateParser::parseISO("2021-03-05"), DateParser::parseISO("2021-06-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-09-06"), notional, DateParser::parseISO("2021-06-07"), DateParser::parseISO("2021-09-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-12-06"), notional, DateParser::parseISO("2021-09-06"), DateParser::parseISO("2021-12-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-03-07"), notional, DateParser::parseISO("2021-12-06"), DateParser::parseISO("2022-03-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-06-07"), notional, DateParser::parseISO("2022-03-07"), DateParser::parseISO("2022-06-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-09-05"), notional, DateParser::parseISO("2022-06-07"), DateParser::parseISO("2022-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-12-05"), notional, DateParser::parseISO("2022-09-05"), DateParser::parseISO("2022-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-03-06"), notional, DateParser::parseISO("2022-12-05"), DateParser::parseISO("2023-03-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-06-05"), notional, DateParser::parseISO("2023-03-06"), DateParser::parseISO("2023-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-09-05"), notional, DateParser::parseISO("2023-06-05"), DateParser::parseISO("2023-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-12-05"), notional, DateParser::parseISO("2023-09-05"), DateParser::parseISO("2023-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-03-05"), notional, DateParser::parseISO("2023-12-05"), DateParser::parseISO("2024-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-06-05"), notional, DateParser::parseISO("2024-03-05"), DateParser::parseISO("2024-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-09-05"), notional, DateParser::parseISO("2024-06-05"), DateParser::parseISO("2024-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-12-05"), notional, DateParser::parseISO("2024-09-05"), DateParser::parseISO("2024-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-03-05"), notional, DateParser::parseISO("2024-12-05"), DateParser::parseISO("2025-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-06-05"), notional, DateParser::parseISO("2025-03-05"), DateParser::parseISO("2025-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-09-05"), notional, DateParser::parseISO("2025-06-05"), DateParser::parseISO("2025-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-12-05"), notional, DateParser::parseISO("2025-09-05"), DateParser::parseISO("2025-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-03-05"), notional, DateParser::parseISO("2025-12-05"), DateParser::parseISO("2026-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-06-05"), notional, DateParser::parseISO("2026-03-05"), DateParser::parseISO("2026-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-09-07"), notional, DateParser::parseISO("2026-06-05"), DateParser::parseISO("2026-09-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-12-07"), notional, DateParser::parseISO("2026-09-07"), DateParser::parseISO("2026-12-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-03-05"), notional, DateParser::parseISO("2026-12-07"), DateParser::parseISO("2027-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-06-07"), notional, DateParser::parseISO("2027-03-05"), DateParser::parseISO("2027-06-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-09-06"), notional, DateParser::parseISO("2027-06-07"), DateParser::parseISO("2027-09-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-12-06"), notional, DateParser::parseISO("2027-09-06"), DateParser::parseISO("2027-12-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-03-06"), notional, DateParser::parseISO("2027-12-06"), DateParser::parseISO("2028-03-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-06-05"), notional, DateParser::parseISO("2028-03-06"), DateParser::parseISO("2028-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-09-05"), notional, DateParser::parseISO("2028-06-05"), DateParser::parseISO("2028-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-12-05"), notional, DateParser::parseISO("2028-09-05"), DateParser::parseISO("2028-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-03-05"), notional, DateParser::parseISO("2028-12-05"), DateParser::parseISO("2029-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-06-05"), notional, DateParser::parseISO("2029-03-05"), DateParser::parseISO("2029-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-09-05"), notional, DateParser::parseISO("2029-06-05"), DateParser::parseISO("2029-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-12-05"), notional, DateParser::parseISO("2029-09-05"), DateParser::parseISO("2029-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-03-05"), notional, DateParser::parseISO("2029-12-05"), DateParser::parseISO("2030-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-06-05"), notional, DateParser::parseISO("2030-03-05"), DateParser::parseISO("2030-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-09-05"), notional, DateParser::parseISO("2030-06-05"), DateParser::parseISO("2030-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-12-05"), notional, DateParser::parseISO("2030-09-05"), DateParser::parseISO("2030-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-03-05"), notional, DateParser::parseISO("2030-12-05"), DateParser::parseISO("2031-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-06-05"), notional, DateParser::parseISO("2031-03-05"), DateParser::parseISO("2031-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-09-05"), notional, DateParser::parseISO("2031-06-05"), DateParser::parseISO("2031-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-12-05"), notional, DateParser::parseISO("2031-09-05"), DateParser::parseISO("2031-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-03-05"), notional, DateParser::parseISO("2031-12-05"), DateParser::parseISO("2032-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-06-07"), notional, DateParser::parseISO("2032-03-05"), DateParser::parseISO("2032-06-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-09-06"), notional, DateParser::parseISO("2032-06-07"), DateParser::parseISO("2032-09-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-12-06"), notional, DateParser::parseISO("2032-09-06"), DateParser::parseISO("2032-12-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-03-07"), notional, DateParser::parseISO("2032-12-06"), DateParser::parseISO("2033-03-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-06-07"), notional, DateParser::parseISO("2033-03-07"), DateParser::parseISO("2033-06-07"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-09-05"), notional, DateParser::parseISO("2033-06-07"), DateParser::parseISO("2033-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-12-05"), notional, DateParser::parseISO("2033-09-05"), DateParser::parseISO("2033-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-03-06"), notional, DateParser::parseISO("2033-12-05"), DateParser::parseISO("2034-03-06"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-06-05"), notional, DateParser::parseISO("2034-03-06"), DateParser::parseISO("2034-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-09-05"), notional, DateParser::parseISO("2034-06-05"), DateParser::parseISO("2034-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-12-05"), notional, DateParser::parseISO("2034-09-05"), DateParser::parseISO("2034-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-03-05"), notional, DateParser::parseISO("2034-12-05"), DateParser::parseISO("2035-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-06-05"), notional, DateParser::parseISO("2035-03-05"), DateParser::parseISO("2035-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-09-05"), notional, DateParser::parseISO("2035-06-05"), DateParser::parseISO("2035-09-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-12-05"), notional, DateParser::parseISO("2035-09-05"), DateParser::parseISO("2035-12-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-03-05"), notional, DateParser::parseISO("2035-12-05"), DateParser::parseISO("2036-03-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-06-05"), notional, DateParser::parseISO("2036-03-05"), DateParser::parseISO("2036-06-05"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-09-05"), notional, DateParser::parseISO("2036-06-05"), DateParser::parseISO("2036-09-05"), fixingDays, index, gearing, spread)));


		std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2021-09-06"),
			DateParser::parseISO("2022-09-05"),
			DateParser::parseISO("2023-09-05"),
			DateParser::parseISO("2024-09-05"),
			DateParser::parseISO("2025-09-05"),
			DateParser::parseISO("2026-09-07"),
			DateParser::parseISO("2027-09-06"),
			DateParser::parseISO("2028-09-05"),
			DateParser::parseISO("2029-09-05"),
			DateParser::parseISO("2030-09-05"),
			DateParser::parseISO("2031-09-05"),
			DateParser::parseISO("2032-09-06"),
			DateParser::parseISO("2033-09-05"),
			DateParser::parseISO("2034-09-05"),
			DateParser::parseISO("2035-09-05"),
			DateParser::parseISO("2036-09-05")
		};

		index->addFixing(DateParser::parseISO("2021-03-04"), 0.0074);

		boost::shared_ptr<LegExerciseOption> option(new LegExerciseOption(optioin_exdates, fixed_payment_dates));

		// 아 그거는 어떻게 하지...? payoff. is zerocallable payoff
		bool isFixedRec = true;

		boost::shared_ptr<Swap> swap(
			new ZeroCouponCallableSwap(fixedLeg, floatingLeg, option, isFixedRec));

		Size timeSteps = Actual365Fixed().yearFraction(evaluationDate, optioin_exdates.back()) * 365;
		Size samples = simulNum();

		boost::shared_ptr<MCZeroCouponCallableSwapEngine> engine(
			new MCZeroCouponCallableSwapEngine(
				process(),
				timeSteps,
				samples,
				1,
				2048));

		swap->setPricingEngine(engine);

		try
		{
			swap->NPV();



		}
		catch (const std::exception& e)
		{
			std::cout << e.what();
		}
	}

	void test_34370()
	{
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-04-30");
		Date evaluationDate = Settings::instance().evaluationDate();
		// get curve 랑 등등 해서 완성했어.
		// 이거는 여기다가 그냥 hard coding 해서 넣는 방식임. 그래서 여러가지 바꿔보는 것...?
		Real notional = 20000000000;

		Leg fixedLeg;
		Real rate = 0.0217;
		std::vector<Date> fixed_payment_dates
		{
			//DateParser::parseISO("2017-08-23"),
			//DateParser::parseISO("2018-08-23"),
			//DateParser::parseISO("2019-08-23"),
			//DateParser::parseISO("2020-08-24"),
			DateParser::parseISO("2022-01-21"),
			DateParser::parseISO("2023-01-23"),
			DateParser::parseISO("2024-01-22"),
			DateParser::parseISO("2025-01-21"),
			DateParser::parseISO("2026-01-21"),
			DateParser::parseISO("2027-01-21"),
			DateParser::parseISO("2028-01-21"),
			DateParser::parseISO("2029-01-22"),
			DateParser::parseISO("2030-01-21"),
			DateParser::parseISO("2031-01-21"),
			DateParser::parseISO("2032-01-21"),
			DateParser::parseISO("2033-01-21"),
			DateParser::parseISO("2034-01-23"),
			DateParser::parseISO("2035-01-22"),
			DateParser::parseISO("2036-01-21"),
			DateParser::parseISO("2037-01-21"),
			DateParser::parseISO("2038-01-21"),
			DateParser::parseISO("2039-01-21"),
			DateParser::parseISO("2040-01-23")
		};

		std::vector<Real> amounts{
			877417800,
			1330457766,
			1793328700,
			2266243933,
			2749421426,
			3243083871,
			3747458791,
			4262778647,
			4789280943,
			5327208340,
			5876808761,
			6438335511,
			7012047391,
			7598208820,
			8197089951,
			8808966803,
			9434121383,
			10072841817,
			10725422484
		};

		for (Size i = 0; i < fixed_payment_dates.size(); i++)
		{
			boost::shared_ptr<SimpleCashFlow> scf(new SimpleCashFlow(amounts[i], fixed_payment_dates[i]));
			fixedLeg.push_back(scf);
		}

		//for (Size i = 0; i < fixed_payment_dates.size(); i++)
		//{
		//	Real amount = notional * rate * (i + 2);
		//	boost::shared_ptr<SimpleCashFlow> scf(new SimpleCashFlow(amount, fixed_payment_dates[i]));
		//	fixedLeg.push_back(scf);
		//}

		Leg floatingLeg;

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new Libor("KRWCD3M", Period(3, TimeUnit::Months), 1,
				KRWCurrency(), SouthKorea(), Actual365Fixed(), yts_h));

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.0;

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-04-21"), notional, DateParser::parseISO("2020-01-21"), DateParser::parseISO("2020-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-07-21"), notional, DateParser::parseISO("2020-04-21"), DateParser::parseISO("2020-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-10-21"), notional, DateParser::parseISO("2020-07-21"), DateParser::parseISO("2020-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-01-21"), notional, DateParser::parseISO("2020-10-21"), DateParser::parseISO("2021-01-21"), fixingDays, index, gearing, spread)));

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-04-21"), notional, DateParser::parseISO("2021-01-21"), DateParser::parseISO("2021-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-07-21"), notional, DateParser::parseISO("2021-04-21"), DateParser::parseISO("2021-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-10-21"), notional, DateParser::parseISO("2021-07-21"), DateParser::parseISO("2021-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-01-21"), notional, DateParser::parseISO("2021-10-21"), DateParser::parseISO("2022-01-21"), fixingDays, index, gearing, spread)));

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-04-21"), notional, DateParser::parseISO("2022-01-21"), DateParser::parseISO("2022-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-07-21"), notional, DateParser::parseISO("2022-04-21"), DateParser::parseISO("2022-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-10-21"), notional, DateParser::parseISO("2022-07-21"), DateParser::parseISO("2022-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-01-23"), notional, DateParser::parseISO("2022-10-21"), DateParser::parseISO("2023-01-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-04-21"), notional, DateParser::parseISO("2023-01-23"), DateParser::parseISO("2023-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-07-21"), notional, DateParser::parseISO("2023-04-21"), DateParser::parseISO("2023-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-10-23"), notional, DateParser::parseISO("2023-07-21"), DateParser::parseISO("2023-10-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-01-22"), notional, DateParser::parseISO("2023-10-23"), DateParser::parseISO("2024-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-04-22"), notional, DateParser::parseISO("2024-01-22"), DateParser::parseISO("2024-04-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-07-22"), notional, DateParser::parseISO("2024-04-22"), DateParser::parseISO("2024-07-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-10-21"), notional, DateParser::parseISO("2024-07-22"), DateParser::parseISO("2024-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-01-21"), notional, DateParser::parseISO("2024-10-21"), DateParser::parseISO("2025-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-04-21"), notional, DateParser::parseISO("2025-01-21"), DateParser::parseISO("2025-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-07-21"), notional, DateParser::parseISO("2025-04-21"), DateParser::parseISO("2025-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-10-21"), notional, DateParser::parseISO("2025-07-21"), DateParser::parseISO("2025-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-01-21"), notional, DateParser::parseISO("2025-10-21"), DateParser::parseISO("2026-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-04-21"), notional, DateParser::parseISO("2026-01-21"), DateParser::parseISO("2026-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-07-21"), notional, DateParser::parseISO("2026-04-21"), DateParser::parseISO("2026-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-10-21"), notional, DateParser::parseISO("2026-07-21"), DateParser::parseISO("2026-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-01-21"), notional, DateParser::parseISO("2026-10-21"), DateParser::parseISO("2027-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-04-21"), notional, DateParser::parseISO("2027-01-21"), DateParser::parseISO("2027-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-07-21"), notional, DateParser::parseISO("2027-04-21"), DateParser::parseISO("2027-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-10-21"), notional, DateParser::parseISO("2027-07-21"), DateParser::parseISO("2027-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-01-21"), notional, DateParser::parseISO("2027-10-21"), DateParser::parseISO("2028-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-04-21"), notional, DateParser::parseISO("2028-01-21"), DateParser::parseISO("2028-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-07-21"), notional, DateParser::parseISO("2028-04-21"), DateParser::parseISO("2028-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-10-23"), notional, DateParser::parseISO("2028-07-21"), DateParser::parseISO("2028-10-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-01-22"), notional, DateParser::parseISO("2028-10-23"), DateParser::parseISO("2029-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-04-23"), notional, DateParser::parseISO("2029-01-22"), DateParser::parseISO("2029-04-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-07-23"), notional, DateParser::parseISO("2029-04-23"), DateParser::parseISO("2029-07-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-10-22"), notional, DateParser::parseISO("2029-07-23"), DateParser::parseISO("2029-10-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-01-21"), notional, DateParser::parseISO("2029-10-22"), DateParser::parseISO("2030-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-04-22"), notional, DateParser::parseISO("2030-01-21"), DateParser::parseISO("2030-04-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-07-22"), notional, DateParser::parseISO("2030-04-22"), DateParser::parseISO("2030-07-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-10-21"), notional, DateParser::parseISO("2030-07-22"), DateParser::parseISO("2030-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-01-21"), notional, DateParser::parseISO("2030-10-21"), DateParser::parseISO("2031-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-04-21"), notional, DateParser::parseISO("2031-01-21"), DateParser::parseISO("2031-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-07-21"), notional, DateParser::parseISO("2031-04-21"), DateParser::parseISO("2031-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-10-21"), notional, DateParser::parseISO("2031-07-21"), DateParser::parseISO("2031-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-01-21"), notional, DateParser::parseISO("2031-10-21"), DateParser::parseISO("2032-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-04-21"), notional, DateParser::parseISO("2032-01-21"), DateParser::parseISO("2032-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-07-21"), notional, DateParser::parseISO("2032-04-21"), DateParser::parseISO("2032-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-10-21"), notional, DateParser::parseISO("2032-07-21"), DateParser::parseISO("2032-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-01-21"), notional, DateParser::parseISO("2032-10-21"), DateParser::parseISO("2033-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-04-21"), notional, DateParser::parseISO("2033-01-21"), DateParser::parseISO("2033-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-07-21"), notional, DateParser::parseISO("2033-04-21"), DateParser::parseISO("2033-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2033-10-21"), notional, DateParser::parseISO("2033-07-21"), DateParser::parseISO("2033-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-01-23"), notional, DateParser::parseISO("2033-10-21"), DateParser::parseISO("2034-01-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-04-21"), notional, DateParser::parseISO("2034-01-23"), DateParser::parseISO("2034-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-07-21"), notional, DateParser::parseISO("2034-04-21"), DateParser::parseISO("2034-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2034-10-23"), notional, DateParser::parseISO("2034-07-21"), DateParser::parseISO("2034-10-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-01-22"), notional, DateParser::parseISO("2034-10-23"), DateParser::parseISO("2035-01-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-04-23"), notional, DateParser::parseISO("2035-01-22"), DateParser::parseISO("2035-04-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-07-23"), notional, DateParser::parseISO("2035-04-23"), DateParser::parseISO("2035-07-23"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2035-10-22"), notional, DateParser::parseISO("2035-07-23"), DateParser::parseISO("2035-10-22"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-01-21"), notional, DateParser::parseISO("2035-10-22"), DateParser::parseISO("2036-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-04-21"), notional, DateParser::parseISO("2036-01-21"), DateParser::parseISO("2036-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-07-21"), notional, DateParser::parseISO("2036-04-21"), DateParser::parseISO("2036-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2036-10-21"), notional, DateParser::parseISO("2036-07-21"), DateParser::parseISO("2036-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2037-01-21"), notional, DateParser::parseISO("2036-10-21"), DateParser::parseISO("2037-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2037-04-21"), notional, DateParser::parseISO("2037-01-21"), DateParser::parseISO("2037-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2037-07-21"), notional, DateParser::parseISO("2037-04-21"), DateParser::parseISO("2037-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2037-10-21"), notional, DateParser::parseISO("2037-07-21"), DateParser::parseISO("2037-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2038-01-21"), notional, DateParser::parseISO("2037-10-21"), DateParser::parseISO("2038-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2038-04-21"), notional, DateParser::parseISO("2038-01-21"), DateParser::parseISO("2038-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2038-07-21"), notional, DateParser::parseISO("2038-04-21"), DateParser::parseISO("2038-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2038-10-21"), notional, DateParser::parseISO("2038-07-21"), DateParser::parseISO("2038-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2039-01-21"), notional, DateParser::parseISO("2038-10-21"), DateParser::parseISO("2039-01-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2039-04-21"), notional, DateParser::parseISO("2039-01-21"), DateParser::parseISO("2039-04-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2039-07-21"), notional, DateParser::parseISO("2039-04-21"), DateParser::parseISO("2039-07-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2039-10-21"), notional, DateParser::parseISO("2039-07-21"), DateParser::parseISO("2039-10-21"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2040-01-23"), notional, DateParser::parseISO("2039-10-21"), DateParser::parseISO("2040-01-23"), fixingDays, index, gearing, spread)));



		std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2022-01-21"),
			DateParser::parseISO("2023-01-23"),
			DateParser::parseISO("2024-01-22"),
			DateParser::parseISO("2025-01-21"),
			DateParser::parseISO("2026-01-21"),
			DateParser::parseISO("2027-01-21"),
			DateParser::parseISO("2028-01-21"),
			DateParser::parseISO("2029-01-22"),
			DateParser::parseISO("2030-01-21"),
			DateParser::parseISO("2031-01-21"),
			DateParser::parseISO("2032-01-21"),
			DateParser::parseISO("2033-01-21"),
			DateParser::parseISO("2034-01-23"),
			DateParser::parseISO("2035-01-22"),
			DateParser::parseISO("2036-01-21"),
			DateParser::parseISO("2037-01-21"),
			DateParser::parseISO("2038-01-21"),
			DateParser::parseISO("2039-01-21"),
			DateParser::parseISO("2040-01-23")
		};

		index->addFixing(DateParser::parseISO("2021-04-20"), 0.0073);

		boost::shared_ptr<LegExerciseOption> option(new LegExerciseOption(optioin_exdates, fixed_payment_dates));

		// 아 그거는 어떻게 하지...? payoff. is zerocallable payoff
		bool isFixedRec = true;

		boost::shared_ptr<Swap> swap(
			new ZeroCouponCallableSwap(fixedLeg, floatingLeg, option, isFixedRec));

		Size timeSteps = Actual365Fixed().yearFraction(evaluationDate, optioin_exdates.back()) * 365;
		Size samples = simulNum();

		boost::shared_ptr<MCZeroCouponCallableSwapEngine> engine(
			new MCZeroCouponCallableSwapEngine(
				process(),
				timeSteps,
				samples,
				1,
				2048));

		swap->setPricingEngine(engine);

		try
		{
			swap->NPV();
			// a360    : -2,470,225,814
			// calypso : -2,563,845,632
			//


		}
		catch (const std::exception& e)
		{
			std::cout << e.what();
		}
	}

	void test_from_input()
	{
		std::vector<Real> alpha_times{ 100.0 };
		PiecewiseConstantParameter alpha_para(alpha_times);
		alpha_para.setParam(0, 0.05);

		std::vector<Real> sigma_times{ 0.083333,0.25,0.5,0.75,1,1.5,2,3,4,5,7,10 };
		std::vector<Real> sigma_values
		{
			0.006188,
			0.00580129905905104,
			0.00566797309450212,
			0.00539632208082505,
			0.00509562518244817,
			0.00490434858059661,
			0.00438993086505927,
			0.00461860790715125,
			0.00428171893052311,
			0.00393527394217989,
			0.004635,
			0.004635
		};

		LinearInterpolationParameter sigma_para(sigma_times, sigma_values);

		for (Size i = 0; i < sigma_times.size(); i++)
		{
			std::cout << sigma_times[i] << " , " << sigma_para(sigma_times[i]) << std::endl;
		}

		std::cout << "-------------------------" << std::endl;

		for (Size i = 0; i < 100; i++)
		{
			std::cout << i * 0.1 << " , " << sigma_para(i*0.1) << std::endl;
		}

		std::cout << sigma_times.back() - 0.1 << " , " << sigma_para(sigma_times.back() - 0.1) << std::endl;
		std::cout << sigma_times.back() << " , " << sigma_para(sigma_times.back()) << std::endl;
		std::cout << sigma_times.back() + 0.1 << " , " << sigma_para(sigma_times.back() + 0.1) << std::endl;

		// 여기서 한꺼번에 집어 넣는 작업

		// read from json

	}

	void get_sample_input()
	{
		// 테스트는...?!

		// save to json file
	}
}


// 그냥 여기서 다 처리함
//void calculate3() const
//{
//	std::vector < boost::function1<Real, Real> v_ =
//		LsmBasisSystem::pathBasisSystem(polynomOrder, polynomType));

//	// 이거는 어떻게 정하지....?
//	// call schedule이 coupon schedule의 subset 이어야함.
//	Leg fixedLeg;
//	Leg floatingLeg;

//	std::vector<boost::shared_ptr<LegExerciseOption>> options;

//	Date evaluationDate = Settings::instance().evaluationDate();
//	Actual365Fixed daycounter;
//	// 다 필요없고 그냥 월단위나 일단위로 뽑을가...? 월로 해야함.
//	// 월은 effective 기준으로


//	// calibration phase
//	// state별 가격을 전부 넣어야함.

//	// calibration phase
//	// 돌면서 coeff가 있어야함
//	// 그 coeff로 -> 가격을 추정함

//	const Size len = options.size();
//	const Size n = this->samples_;
//	Array prices(n), exercise(n);

//	std::vector<Time> excerciseTimes;
//	for (Size i = 0; i < len; i++)
//	{
//		Time t = daycounter.yearFraction(evaluationDate, exDate);
//		Date exDate = options[i]->lastDate();

//		if (exDate < evaluationDate)
//			continue;

//		excerciseTimes.push_back(t);
//	}

//	// [t_len][simul]
//	std::vector<std::vector<Real>> x_items(len); // []
//	std::vector<std::vector<Real>> y_items(len);
//	std::vector<std::vector<Real>> discounts(len);

//	for (Size j = 0; j < n; j++)
//	{
//		Path path;

//		// | ----- | ----- | ----- | ----- |
//		// | ----- | ----- | ----- | ----- |
//		// | ----- | ----- | ----- | ----- |
//		// | ----- | ----- | ----- | ----- |
//		// std::vector<Real> discounts = this->discounts(path, excerciseTimes);

//		// accumulated
//		std::vector<Real> float_amt = this->legFltAmt(fixedLeg, path, excerciseTimes);
//		std::vector<Real> fixed_amt = this->legFixAmt(floatingLeg, path, excerciseTimes);
//		std::vector<Real> states = this->states(path, excerciseTimes);

//		Real npv = fixed_amt - float_amt;

//		std::vector<Real>& x = x_items[j];
//		std::vector<Real>& y = y_items[j];

//		// roll back step
//		for (Size i = 0; i < len; ++i) {
//			exercise[j] = float_amt[i] - fixed_amt[i];

//
//			x.push_back(states[i]);
//			y.push_back(dF_[i] * prices[i]);

//			//if (exercise[j] > 0.0) {
//			//}
//		}
//	}

//	for (Size i = 0; i < len; ++i)
//	{
//		if (v_.size() <= x.size()) {
//			coeff_[i] = GeneralLinearLeastSquares(x, y, v_).coefficients();
//		}
//		else {
//			// if number of itm paths is smaller then the number of
//			// calibration functions then early exercise if exerciseValue > 0
//			coeff_[i] = Array(v_.size(), 0.0);
//		}
//	}

//	for (Size j = 0, k = 0; j < n; ++j)
//	{
//		std::vector<Real>& x = x_items[j];
//		std::vector<Real>& y = y_items[j];

//		for (Size i = 0; i < len; ++i)
//		{
//			prices[j] *= discounts[i][j]; // 한번 땡김

//			Real continuationValue = 0.0;

//			for (Size l = 0; l < v_.size(); ++l) {
//				continuationValue += coeff_[i][l] * v_[l](x[k]);
//			}
//			if (continuationValue < exercise[j]) {
//				prices[j] = exercise[j];
//			}
//			++k;

//			//if (exercise[j] > 0.0) {
//			//}
//		}
//	}

//	// backward
//	// | ----- | ----- | ----- | ----- |
//	//	                       r
//	//                                 v[i]
//	//                         v[i-1] = v[i] * df(accumulated shortrates)
//	//
//	for (Size i = len - 1; i > 0; --i)
//	{

//	}

//	// ex 부분이있고, 각각의 payoff부분이 있어.

//	// pricing phase

//	// simulation 별로 구함.

//	// 보니까 순서가
//	// 처음에 calibration 을 하고

//}


/*
// t 기준의 amount // float 부리 해야함.
std::vector<Real> legFltAmt(const Leg& leg, const Path& path, std::vector<Date> excerciseDates) const
{
	Date evaluationDate = Settings::instance().evaluationDate();

	Size len = excerciseDates.size();

	std::vector<Size> ex_t_pos;
	std::vector<Real> amounts;
	std::vector<Real> rates;

	Size legNum = leg.size();
	Size startLegCount = 0;

	for (Size i = 0; i < legNum; i++)
	{
		startLegCount = i;
		Date d = leg[i]->date();
		if (evaluationDate < d) // 당일꺼는 안넣음
			break;
	}

	for (Size i = 0; i < len; i++)
	{
		Date excerciseDate = excerciseDates[i];
		if (excerciseDate <= evaluationDate) // 당일꺼는 안넣음
			continue;

		Real amount = 0.0;

		for (Size j = startLegCount; j < legNum; j++)
		{
			boost::shared_ptr<FloatingRateCoupon> fltCpn
				= boost::dynamic_pointer_cast<FloatingRateCoupon>(leg[j]);
			Date paymentDate = fltCpn->date();
			if (excerciseDate < paymentDate)
				break;

			DayCounter dc = fltCpn->dayCounter();

			Real t = dc.yearFraction(evaluationDate, paymentDate);
			Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

			if (reset_t < 0)
			{
				amount += fltCpn->amount();
			}
			else
			{
				Real notional = fltCpn->notional();
				Size pos = path.timeGrid().closestIndex(reset_t);
				Real rate = this->rate(path, reset_t);
				Real accrualPeriod = fltCpn->accrualPeriod();
				Real compound = 1.0 / this->hw_model_->discountBond(reset_t, reset_t + accrualPeriod, rate);
				Real float_rate = std::pow(compound, 1.0 / accrualPeriod)- 1.0;

				Real coupon = notional * float_rate * accrualPeriod; // 이게 그날 amount;
				Real toExerciseDate_accrualPeriod = dc.yearFraction(paymentDate, excerciseDate);
				Real coupon_at_exercise = coupon * this->compound(path, paymentDate, excerciseDate, dc); //  다음으로 보냄
				amount += coupon_at_exercise;

				// amount *= exp(this->rate(path, t) * accrualPeriod); //  다음으로 보냄
			}

		}

		amounts.push_back(amount);

	}

	return amounts;
}

// t 기준의 amount
std::vector<Real> legFixAmt(const Leg& leg, const Path& path, std::vector<Date> excerciseDates) const
{
	Date evaluationDate = Settings::instance().evaluationDate();
	Size len = excerciseDates.size();

	std::vector<Real> amounts;

	Size legNum = leg.size();
	Size startLegCount = 0;
	Real amount = 0.0;

	QL_REQUIRE(legNum == len, "fixedLeg size must be equal to excerciseDates size");

	for (Size i = 0; i < legNum; i++)
	{
		startLegCount = i;
		if (evaluationDate < leg[i]->date()) // 당일꺼는 안넣음
			break;
	}

	for (Size i = startLegCount; i < legNum; i++)
	{
		boost::shared_ptr<SimpleCashFlow> simpleCpn
			= boost::dynamic_pointer_cast<SimpleCashFlow>(leg[i]);

		amount = simpleCpn->amount();
		amounts.push_back(amount);
	}

	return amounts;
}



std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2021-07-28"),
			DateParser::parseISO("2022-07-27"),
			DateParser::parseISO("2023-07-27"),
			DateParser::parseISO("2024-07-29"),
			DateParser::parseISO("2025-07-29"),
			DateParser::parseISO("2026-07-29"),
			DateParser::parseISO("2027-07-28"),
			DateParser::parseISO("2028-07-27"),
			DateParser::parseISO("2029-07-27"),
			DateParser::parseISO("2030-07-29"),
			DateParser::parseISO("2031-07-29"),
			DateParser::parseISO("2032-07-28"),
			DateParser::parseISO("2033-07-27"),
			DateParser::parseISO("2034-07-27"),
			DateParser::parseISO("2035-07-27"),
			DateParser::parseISO("2036-08-23") // 만기일
		};

*/