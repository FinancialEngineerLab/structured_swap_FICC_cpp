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
#include <ql/processes/stochasticprocessarray.hpp>
#include <ql/processes/g2extprocess.hpp>
#include <ql/models/shortrate/twofactormodels/g2.hpp>
#include <te/instruments/commons.hpp>
#include <sg/process/model/processmodel.hpp>
#include <qle/factories/xen_class_factory.hpp>
#include <sg/calculations/additionalcalc/ratescalc.hpp>

using namespace QuantLib;
using namespace scenariogenerator;


namespace TestEngineDev
{
	TimeGrid build_timeGrid(const Date& refDate, Size maxYear,
		Size nPerYear,
		const std::string& frequency,
		Size frequency_month, Size frequency_day)
	{
		Time maxYearTime = (double)maxYear;

		QL_REQUIRE(nPerYear > 0, "nPerYear must be positive");
		QL_REQUIRE(nPerYear <= 365, "maximum nPerYear is 365");

		QL_REQUIRE(maxYear <= 200, "maximum maxYear is 200 ");
		QL_REQUIRE(maxYear > 0, "maxYear must be positive");

		QL_REQUIRE(1 <= frequency_month && frequency_month <= 12, "month must between 1 to 12");
		QL_REQUIRE(1 <= frequency_day && frequency_day <= 31, "day must between 1 to 31");

		NullCalendar calendar;
		ActualActual daycounter;

		Size num = static_cast<Size>(maxYear * nPerYear);
		//TimeGrid timeGrid = TimeGrid(maxYearTime, num);

		Date endDate = calendar.advance(refDate, Period(maxYear, TimeUnit::Years));
		Date roopDate = refDate;
		std::vector<Time> times;
		times.push_back(0.0);
		std::vector<Date> timegrid_dates;
		timegrid_dates.clear();
		timegrid_dates.push_back(refDate);

		std::string freq_upp = boost::to_upper_copy(frequency);

		Size pos = freq_upp.find(".");

		if (pos > 0)
			freq_upp = freq_upp.substr(pos + 1);

		bool endOfMonth = false;

		if (freq_upp == "CUSTOM") // date first
		{
			for (Size i = 0; i < maxYear; i++)
			{
				Date next_year_date = calendar.advance(refDate, Period(i, TimeUnit::Years));
				Year next_year = next_year_date.year();
				if (refDate.month() >= 3) // 
					next_year += 1;
				Real actday = Date::isLeap(next_year) ? 366.0 : 365.0;
				Real d_day = actday / (Real)nPerYear;

				for (Size j = 0; j < nPerYear; j++)
				{
					//Size day = static_cast<unsigned int>(d_day)*(j + 1);
					Size day = static_cast<unsigned int>(d_day * (j + 1));
					Date date = calendar.advance(next_year_date, Period(day, TimeUnit::Days));
					timegrid_dates.push_back(date);
					times.push_back(daycounter.yearFraction(refDate, date));
				}
			}

			times.pop_back();
			timegrid_dates.pop_back();

			times.push_back(maxYear * 1.0);
			timegrid_dates.push_back(endDate);

		}
		else if (freq_upp == "MENUAL")
		{

		}
		else
		{
			//Size nPerYear = 1; //Size step = 1;
			Period peoriod(Frequency::Daily);

			if (freq_upp == "DAY" || freq_upp == "DAILY") {
				peoriod = Period(Frequency::Daily);
				times.reserve((maxYear + 1) * 365);
				timegrid_dates.reserve((maxYear + 1) * 365);
			}
			// not use
			else if (freq_upp == "WEEK" || freq_upp == "WEEKLY") { peoriod = Period(Frequency::Weekly); }
			else if (freq_upp == "MONTH" || freq_upp == "QUARTER" || freq_upp == "SEMIANNUAL" ||
				freq_upp == "MONTHLY" || freq_upp == "QUARTERLY" || freq_upp == "SEMIANNUALLY")
			{
				if (freq_upp == "MONTH" || freq_upp == "MONTHLY") { peoriod = Period(Frequency::Monthly); }
				else if (freq_upp == "QUARTER" || freq_upp == "QUARTERLY") { peoriod = Period(Frequency::Quarterly); }
				else if (freq_upp == "SEMIANNUAL" || freq_upp == "SEMIANNUALLY") { peoriod = Period(Frequency::Semiannual); }
				else {}

				bool isLeap = Date::isLeap(roopDate.year());
				while (Date::monthLength(roopDate.month(), isLeap) < frequency_day)
				{
					Date roopDateEndOfMonth = Date::endOfMonth(roopDate);
					timegrid_dates.push_back(roopDateEndOfMonth);
					times.push_back(daycounter.yearFraction(refDate, roopDateEndOfMonth));
					roopDate = calendar.advance(roopDate, peoriod);
					isLeap = Date::isLeap(roopDate.year());
				}

				roopDate = Date(frequency_day, roopDate.month(), roopDate.year());
				endDate = calendar.advance(roopDate, Period(maxYear, TimeUnit::Years));

				if (refDate < roopDate)
				{
					timegrid_dates.push_back(roopDate);
					times.push_back(daycounter.yearFraction(refDate, roopDate));
				}

			}
			else if (freq_upp == "ANNUAL" || freq_upp == "ANNUALLY") {
				peoriod = Period(Frequency::Annual);

				roopDate = Date(frequency_day, Month(frequency_month), roopDate.year());

				if (refDate < roopDate)
				{
					timegrid_dates.push_back(roopDate);
					times.push_back(daycounter.yearFraction(refDate, roopDate));
				}

			}
			else if (freq_upp == "FIRSTOFMONTH") {
				peoriod = Period(Frequency::Monthly);
				roopDate = Date(1, roopDate.month(), roopDate.year());
			}
			else if (freq_upp == "FIRSTOFANNUAL") {
				peoriod = Period(Frequency::Annual);
				roopDate = Date(1, Month(1), roopDate.year() + 1);
				timegrid_dates.push_back(roopDate);
				times.push_back(daycounter.yearFraction(refDate, roopDate));
			}
			else if (freq_upp == "ENDOFMONTH") {
				peoriod = Period(Frequency::Monthly);
				endOfMonth = true;
				Date thisendofMonth = Date::endOfMonth(roopDate);
				if (thisendofMonth != roopDate)
				{
					timegrid_dates.push_back(thisendofMonth);
					times.push_back(daycounter.yearFraction(refDate, thisendofMonth));
				}
			}
			else if (freq_upp == "ENDOFANNUAL") {
				endOfMonth = true;
				peoriod = Period(Frequency::Annual);
				Date thisendofAnnual = Date::endOfMonth(Date(1, Month(12), roopDate.year()));

				if (thisendofAnnual != roopDate)
				{
					timegrid_dates.push_back(thisendofAnnual);
					times.push_back(daycounter.yearFraction(refDate, thisendofAnnual));
					roopDate = thisendofAnnual;
				}
			}
			else {
				QL_FAIL("unknown timegrid frequency : " << freq_upp << "\n"
					<< "day, week, month, semiannual, annual, firstofmonth, firstofquarter, firstofsemiannual, firstofannual, endofmonth, endofquarter, endofsemiannual, endofannual");
			}

			int roopCount = 1;
			Date startDate = roopDate;

			//std::cout << "timedate grid building...  ";

			while (roopDate < endDate)
			{
				//roopDate = calendar.advance(roopDate, peoriod);
				//roopDate = calendar.advance(startDate, roopCount * peoriod);
				roopDate = startDate + roopCount * peoriod;
				if (endOfMonth)
					roopDate = Date::endOfMonth(roopDate);
				timegrid_dates.push_back(roopDate);
				times.push_back(daycounter.yearFraction(refDate, roopDate));
				roopCount += 1;
				//std::cout << ".";
				//if (roopCount % 100 == 0)
				//	std::cout << std::endl;
			}

			//std::cout << "done." << std::endl;

			//if (roopDate != endDate)
			if (roopDate < endDate)
				timegrid_dates.push_back(endDate);

		}

		return TimeGrid(times.begin(), times.end(), timegrid_dates);
	}

	// CMS Spread
	class StructuredSwap : public Swap
	{
	public:
		class arguments : public virtual PricingEngine::arguments {
		public:
			Leg payLeg;
			Leg recLeg;
			std::vector<Date> excerciseDates;
			void validate() const {}
		};

		class results : public Instrument::results {
		public:
			void reset() {}
		};

		class engine : public GenericEngine<StructuredSwap::arguments,
			StructuredSwap::results> {};
		//ZeroCouponCallableSwap(
		//	std::vector<Leg> legs,
		//	const std::vector<bool>& payer,
		//	const std::vector<boost::shared_ptr<LegExerciseOption>> options)
		//: Swap(legs, payer), options_(options)
		//{
		//}

		// zero callable 이니까 fixed 기준으로 함.
		// options 에 들어있는 ex date는 floatingLeg의 subset이어야함?
		StructuredSwap(
			Leg payLeg,
			Leg recLeg,
			const boost::shared_ptr<LegExerciseOption>& option)
			: Swap(payLeg, recLeg), option_(option)
		{
			//boost::shared_ptr<IborCouponPricer> pricer(new BlackIborCouponPricer);
			//setCouponPricer(floatingLeg, pricer);
		}

		void setupArguments(PricingEngine::arguments* args) const
		{
			StructuredSwap::arguments* arguments = dynamic_cast<StructuredSwap::arguments*>(args);
			QL_REQUIRE(arguments != 0, "wrong argument type");

			arguments->excerciseDates = this->option_->dates();
			arguments->payLeg = this->leg(0);
			arguments->recLeg = this->leg(1);
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

	class ScenarioPath
	{
	public:
		// ScenarioPath(const MultiPath& multiPath, std::vector<Array> factors_arr) : multiPath(multiPath), factors_arr(factors_arr)
		ScenarioPath(const std::vector<std::vector<Real>>& data, const TimeGrid& timeGrid)
			: data_(data), timeGrid_(timeGrid)
		{
			this->assetNum_ = data.size();
			this->current_underlyings = std::vector<Real>(this->assetNum_);
		}

		void set_current_underlyings(Size pos) const
		{
			for (Size i = 0; i < this->assetNum_; i++)
				this->current_underlyings[i] = data_[i][pos];
		}

	public:
		// const MultiPath& multiPath;
		//std::vector<Path> data; // value of indexes

		// mutable std::vector<Array> factors_arr; // 각 모델별 factor들임.
		mutable std::vector<Real> current_underlyings; // 각 underlying 별 값들. 
		mutable Time t;
		mutable Size pos;

	private:
		std::vector<std::vector<Real>> data_;; // value of indexes
		const TimeGrid& timeGrid_;
		Size assetNum_;
	};

	class PayoffMC
	{
	public:
		virtual ~PayoffMC() {}
		virtual std::string name() const = 0;
		virtual Real calc_payoff(const ScenarioPath& scenPath) const = 0;
	};

	// 우선 cms에 대해서 함
	// 여기 안에다가 factor 랑 설정을 함.
	// 그러고 무조건 multiPath를 받음.
	class ConditionMC
	{
	public:
		virtual ~ConditionMC() {}
		virtual bool check(const ScenarioPath& scenPath) const = 0;

	protected:

	};

	class AndConditionMC : public ConditionMC
	{
	public:
		AndConditionMC(const std::vector<boost::shared_ptr<ConditionMC>>& conditions)
			: conditions_(conditions) { }

		bool check(const ScenarioPath& scenPath) const
		{
			for (Size i = 0; i < this->conditions_.size(); i++)
			{
				if (!this->conditions_[i]->check(scenPath))
					return false;
			}
			return true;
		}

	protected:
		std::vector<boost::shared_ptr<ConditionMC>> conditions_;
	};

	class RangeConditionMC : public ConditionMC
	{
	public:
		RangeConditionMC(const boost::shared_ptr<PayoffMC>& po,
			Real a, Real b)
			: a_(a), b_(b) {}
		// bool check(const ScenarioPath& scenPath) const
		bool check(const ScenarioPath& scenPath) const
		{
			Real v = po_->calc_payoff(scenPath);
			if (a_ <= v && v <= b_)
				return true;
			return false;
		}

	protected:
		Real a_;
		Real b_;
		boost::shared_ptr<PayoffMC> po_;
	};


	class LinearPayoffMC : public PayoffMC
	{
	public:
		LinearPayoffMC(const boost::shared_ptr<PayoffMC>& po,
			Real multiple, Real spread) {}
		virtual Real calc_payoff(const ScenarioPath& scenPath) const = 0;
	};

#pragma region Operator

	// unary
	boost::shared_ptr<PayoffMC> operator+(const boost::shared_ptr<PayoffMC>& po);
	boost::shared_ptr<PayoffMC> operator-(const boost::shared_ptr<PayoffMC>& po);

	// binary
	boost::shared_ptr<PayoffMC> operator+(const boost::shared_ptr<PayoffMC>&, const boost::shared_ptr<PayoffMC>&);
	boost::shared_ptr<PayoffMC> operator-(const boost::shared_ptr<PayoffMC>&, const boost::shared_ptr<PayoffMC>&);
	boost::shared_ptr<PayoffMC> operator*(const boost::shared_ptr<PayoffMC>&, const boost::shared_ptr<PayoffMC>&);
	boost::shared_ptr<PayoffMC> operator/(const boost::shared_ptr<PayoffMC>&, const boost::shared_ptr<PayoffMC>&);

	// binary const
	boost::shared_ptr<PayoffMC> operator+(const boost::shared_ptr<PayoffMC>& po, Real v);
	boost::shared_ptr<PayoffMC> operator-(const boost::shared_ptr<PayoffMC>& po, Real v);
	boost::shared_ptr<PayoffMC> operator*(const boost::shared_ptr<PayoffMC>& po, Real v);
	boost::shared_ptr<PayoffMC> operator/(const boost::shared_ptr<PayoffMC>& po, Real v);

	// binary const reverse
	boost::shared_ptr<PayoffMC> operator+(Real v, const boost::shared_ptr<PayoffMC>& po);
	boost::shared_ptr<PayoffMC> operator-(Real v, const boost::shared_ptr<PayoffMC>& po);
	boost::shared_ptr<PayoffMC> operator*(Real v, const boost::shared_ptr<PayoffMC>& po);
	boost::shared_ptr<PayoffMC> operator/(Real v, const boost::shared_ptr<PayoffMC>& po);

#pragma endregion

	class ConditionPayoffMC : public PayoffMC
	{
	public:
		ConditionPayoffMC(const boost::shared_ptr<ConditionMC>& conditionMC,
			const boost::shared_ptr<PayoffMC>& payoffMCTrue,
			const boost::shared_ptr<PayoffMC>& payoffMCFalse)
			: conditionMC_(conditionMC), payoffMCTrue_(payoffMCTrue), payoffMCFalse_(payoffMCFalse)
		{
		}

		Real calc_payoff(const ScenarioPath& scenPath) const
		{
			if (this->conditionMC_->check(scenPath))
				return this->payoffMCTrue_->calc_payoff(scenPath);
			else
				return this->payoffMCFalse_->calc_payoff(scenPath);
		}

	protected:
		boost::shared_ptr<ConditionMC> conditionMC_;
		boost::shared_ptr<PayoffMC> payoffMCTrue_;
		boost::shared_ptr<PayoffMC> payoffMCFalse_;
	};


	// indexes ------------------------------------------------------------------------------------------------------------------------
	class IndexMC : public Index, public PayoffMC
	{
	public:
		IndexMC(const boost::shared_ptr<Index>& index,
			const boost::shared_ptr<ProcessModel>& model,
			const TimeGrid& timeGrid)
			: index_(index), model_(model)
		{
		}

		std::string name() const { return this->index_->name(); }
		Calendar fixingCalendar() const { return this->index_->fixingCalendar(); }
		bool isValidFixingDate(const Date& fixingDate) const { return this->index_->isValidFixingDate(fixingDate); }
		Rate fixing(const Date& fixingDate, bool forecastTodaysFixing = false) const {
			return this->index_->fixing(fixingDate, forecastTodaysFixing);
		}
		boost::shared_ptr<ProcessModel> model() const { return this->model_; }

		// model 이 두개가 들어오는 경우가 있나...? -> 없는 듯 . 하나에서 두개를 뽑는 건 있어도...
		//virtual Real calc_index(Real t, const Array& factors, const boost::shared_ptr<AffineModel>& model) const = 0;
		Real calc_payoff(const ScenarioPath& scenPath) const { return scenPath.current_underlyings[this->location_]; }
		void setLocation(Size location) const { this->location_ = location; } // model location in scenPath
		Size getLocation() const { return this->location_; }

		virtual boost::shared_ptr<AdditionalCalc> get_additionalCalc() const = 0;

	protected:
		boost::shared_ptr<Index> index_;
		boost::shared_ptr<ProcessModel> model_;
		TimeGrid timeGrid_;
		mutable Size location_; // model location in scenPath

	};

	class SwapIndexMC : public IndexMC
	{
	public:
		SwapIndexMC(const boost::shared_ptr<SwapIndex>& swapIndex,
			const boost::shared_ptr<ProcessModel>& model,
			const TimeGrid& timeGrid)
			: IndexMC(swapIndex, model, timeGrid)
		{
			Date evaluationDate = Settings::instance().evaluationDate();
			Date effectiveDate = swapIndex->fixingCalendar().advance(
				evaluationDate, Period(swapIndex->fixingDays(), TimeUnit::Days), swapIndex->fixedLegConvention());
			auto swap = swapIndex->underlyingSwap(effectiveDate);
			auto dates = swap->fixedSchedule().dates();
			auto dc = swapIndex->dayCounter();

			auto fixedLeg = swap->fixedLeg();

			for (Size i = 0; i < fixedLeg.size(); i++)
			{
				this->coupon_frac_arr_.push_back(dc.yearFraction(evaluationDate, fixedLeg[i]->date()));
			}

			//// 처음꺼 제외
			//for (Size i = 1; i < dates.size(); i++)
			//{
			//	this->coupon_frac_arr_.push_back(dc.yearFraction(evaluationDate, dates[i]));
			//}
		}


		//Real calc_index(Real t, const Array& factors, const boost::shared_ptr<AffineModel>& model) const
		//{
		//	Real df_T = model->discountBond(t, t + this->coupon_frac_arr_.back(), factors);
		//	Real sum_df = 0.0;

		//	for (Size i = 0; i < this->coupon_frac_arr_.size(); i++)
		//		sum_df += model->discountBond(t, t + this->coupon_frac_arr_[i], factors);

		//	return (1.0 - df_T) / sum_df;
		//}

		boost::shared_ptr<AdditionalCalc> get_additionalCalc() const
		{
			//
			boost::shared_ptr<AdditionalCalc> cms_calc;
			return cms_calc;
		}


		Real calc_payoff(const ScenarioPath& scenPath) const
		{
			return scenPath.current_underlyings[this->location_];
		}


	protected:
		// boost::shared_ptr<SwapIndex> swapIndex_;
		std::vector<Real> coupon_frac_arr_;
	};

	class IborIndexMC : public IndexMC
	{
	public:
		IborIndexMC(const boost::shared_ptr<IborIndex>& iborIndex,
			const boost::shared_ptr<ProcessModel>& model,
			const TimeGrid& timeGrid)
			: IndexMC(iborIndex, model, timeGrid)
		{
			// 이방식은 초기에 정해지는 거라
			// 중간중간 나오는 yearfrac으로 해야하기 때문에
			// 그냥 0.25로...? -> 근데 이거는 swap도 같은 이슈인데...?
			// 

			// PeriodParser::parse_yearfrac(iborIndex->tenor());
			this->tenor_frac_ = 0.25;
		}

		Real calc_index(Real t, const Array& factors, const boost::shared_ptr<AffineModel>& model) const
		{
			Real compound = 1.0 / model->discountBond(t, t + this->tenor_frac_, factors);
			Real float_rate = (compound - 1.0) / this->tenor_frac_;

			return float_rate;
		}

		boost::shared_ptr<AdditionalCalc> get_additionalCalc() const
		{
			boost::shared_ptr<AdditionalCalc> ibor_calc;
			return ibor_calc;
		}

		Real calc_payoff(const ScenarioPath& scenPath) const
		{

		}

	protected:
		Real tenor_frac_;
	};


	// coupons ------------------------------------------------------------------------------------------------------------------------
	class CouponMC : public CashFlow
	{
	public:

		virtual DayCounter dayCounter() const = 0;
		virtual Real calculate_path(const ScenarioPath& scenPath) const = 0;
		virtual std::vector<std::string> index_names() const = 0;
	};


	// 모가 문제지 이거...? ㅡ.ㅡ;;
	// 사용할때 mc로 들어오면 좀 이상할거 같아서리...
	// 원래 mc돌릴때 structured leg 쪽은 그렇다 쳐
	// 이거를 원래 floating cpn을 사용하려고 하는 거임 그래서 그냥 박아 넣으려고, 근데
	// 결국은 floating cpn 은 model 이 연결 되어야함. 반드시. 그러지 않고는 못쓰기 때문에...
	// 기존의 ql_cpn을 사용하는 거는 따로 해결하자.
	class FloatingRateCouponMC : public CouponMC
	{
	public:

		FloatingRateCouponMC(const boost::shared_ptr<FloatingRateCoupon>& fltCpn) // index 에 IndexMC가 박혀있어야함.
		{
		}

		FloatingRateCouponMC(
			const Date& paymentDate,
			Real nominal,
			const Date& startDate,
			const Date& endDate,
			Natural fixingDays,
			const boost::shared_ptr<IborIndexMC>& indexMC,
			const DayCounter& dayCounter,
			Real gearing = 1.0,
			Spread spread = 0.0,
			bool isInArrears = false)
		{
			QL_REQUIRE(gearing != 0, "Null gearing not allowed");

			this->location_ = this->indexMC_->getLocation();
			this->accrualPeriod_ = dayCounter.yearFraction(startDate, endDate);

			// registerWith(indexMC);
			// registerWith(Settings::instance().evaluationDate());
		}

		DayCounter dayCounter() const { return dayCounter_; }

		Date fixingDate() const {
			// if isInArrears_ fix at the end of period
			Date refDate = isInArrears_ ? accrualEndDate_ : accrualStartDate_;
			return indexMC_->fixingCalendar().advance(refDate,
				-static_cast<Integer>(fixingDays_), Days, Preceding);
		}

		Real calculate_path(const ScenarioPath& scenPath) const
		{
			Real rate = scenPath.current_underlyings[this->location_];
			
			return this->nominal_ * (this->gearing_ * rate + this->spread_) * this->accrualPeriod_;
		}

	protected:
		Date paymentDate_;
		Real nominal_;
		Date accrualStartDate_, accrualEndDate_;
		Natural fixingDays_;
		boost::shared_ptr<IborIndexMC> indexMC_;
		DayCounter dayCounter_;
		Real gearing_;
		Spread spread_;
		bool isInArrears_;
	protected:
		Real accrualPeriod_;
		Size accrualStart_pos_, accrualEnd_pos_;
		Size location_;
	};

	// avg 는 coupon 에서 처리 되어야함
	// ex)  StructuredFormulaAvgCouponMC
	// arrear 이런거도 coupon 에서 처리
	class StructuredFormulaAccrualCouponMC : public CashFlow
	{
	public:
		StructuredFormulaAccrualCouponMC(const Date& paymentDate,
			Real nominal,
			const boost::shared_ptr<PayoffMC>& payoffMC,
			const Date& accrualStartDate,
			const Date& accrualEndDate,
			const DayCounter& dayCounter,
			const TimeGrid& timeGrid)
		{
			this->accrualPeriod_ = dayCounter.yearFraction(accrualStartDate, accrualEndDate);
		}

		DayCounter dayCounter() const { return dayCounter_; }

		//Real calculate_path(const MultiPath& multiPath, const boost::shared_ptr<AffineModel>& model) const
		Real calculate_path(const ScenarioPath& scenPath) const
		{

			// scenPath.factors_arr // -> 얘는 여기서 세팅해 -> 안쪽에서 재계산 때문에 여기서 세팅해서 나감
				// index 는...? 안쪽에서 index가 있을 거임. 그거는 언제 계산이 되는가...
				// 최종적으로 index를 만났을 때...? 음... ㄴㄴ
				// 엔진에서 indexes 를 가지고 있으니까 거기서 underlying 을 죄다 계산함. 그리고 order settting을 미리 해두고 
				// 가져다가 씀. 그래야 중복계산이 안될듯.

			// 여기서 세팅을 해야함. 음... multiPath 에서.. 아... 필요한게 indexes 만 있으면 되는데...? 다른 곳에서 discount...?

			// interpolation...? ㄴㄴ -> coupon 에 있는 날짜는 전부 들어가야함
			// 근데 외부에서 scenario가 들어오는 경우도 있기 때문에 ... 음... 
			scenPath.set_current_underlyings(this->accrualStart_pos_);
			Real calc_coupon_start = this->payoffMC_->calc_payoff(scenPath);

			scenPath.set_current_underlyings(this->accrualEnd_pos_);
			Real calc_coupon_end = this->payoffMC_->calc_payoff(scenPath);

			std::vector<Real> coupons{ calc_coupon_start, calc_coupon_end };
			Real count = 0.0;

			// accural part
			for (Size i = accrualStart_pos_ + 1; i < accrualEnd_pos_; i++)
			{
				scenPath.set_current_underlyings(i);
				Real coupon = this->payoffMC_->calc_payoff(scenPath);
				coupons.push_back(coupon);
			}

			// coupon
			Real total_coupons_count = static_cast<Real>(coupons.size());
			Real results = std::accumulate(coupons.begin(), coupons.end(), 0.0) / total_coupons_count;

			if (isExistAccrued_) {
				Date evaluationDate = Settings::instance().evaluationDate();
				Real weight = this->dayCounter_.yearFraction(evaluationDate, this->accrualEndDate_) / this->dayCounter_.yearFraction(this->accrualStartDate_, this->accrualEndDate_);
				results = this->accruedAmount_ + results * weight; // 이거 합계가 원래 최대 받을 수 있는 쿠폰과 맞아야함.

			}

			return this->nominal_ * results * this->accrualPeriod_;
		}

	protected:
		Date paymentDate_;
		Real nominal_;
		Date accrualStartDate_, accrualEndDate_;
		Size accrualStart_pos_, accrualEnd_pos_;
		DayCounter dayCounter_;
		boost::shared_ptr<PayoffMC> payoffMC_;
		Real accrualPeriod_;

		bool isExistAccrued_;
		Real accruedAmount_;
	};


	// type1 : scenario internal generation 
	class MCStructuredSwapType1Engine : public TestEngineDev::StructuredSwap::engine {
	public:
		MCStructuredSwapType1Engine(
			std::vector<boost::shared_ptr<IndexMC>>& indexes,
			QuantLib::Matrix index_corr,
			const boost::shared_ptr<ProcessModel>& pay_discount_model,
			const boost::shared_ptr<ProcessModel>& rec_discount_model,
			TimeGrid timeGrid,
			Size samples,
			BigNatural seed,
			Size nCalibrationSamples)
			: indexes_(indexes), pay_discount_model_(pay_discount_model), rec_discount_model_(rec_discount_model),
			timeGrid_(timeGrid), samples_(samples), seed_(seed), nCalibrationSamples_(nCalibrationSamples)
		{
			// index 에 사용될 모델 3개
			// index에 해당하는 rate(ibor or swap) 3개
			// disocunt 에 사용될 모델 2개
			// discount addcalc 2개
			// 총 10개...! 중
			
			// 자기꺼 찾아서 underlying에다가 넣고 위치 설정

			for (Size i = 0; i < indexes.size(); i++)
				indexes[i]->setLocation(i);



		}


		Real compound(const ScenarioPath& scenPath, Date from, Date to, DayCounter dc) const
		{
			Date evaluationDate = Settings::instance().evaluationDate();

			Real fraction = dc.yearFraction(from, to);
			Real from_t = dc.yearFraction(evaluationDate, from);
			Real to_t = dc.yearFraction(evaluationDate, to);

			const TimeGrid& timeGrid = this->timeGrid_;

			Size from_pos = timeGrid.closestIndex(from_t, -1);
			Size to_pos = timeGrid.closestIndex(to_t, -1);

			if (from_pos == to_pos)
				return this->rate(scenPath, from_t);

			Real sum = 0.0;
			for (Size i = from_pos; i < to_pos; i++)
			{
				Real r = this->phi_(this->timeGrid_.at(i)) + multiPath[0][i] + multiPath[1][i];
				sum += r;
			}

			sum = sum / (to_pos - from_pos);

			return exp(sum * fraction);

		}

		Real discount(const ScenarioPath& scenPath, Date from, Date to, DayCounter dc) const
		{
			return 1.0 / compound(multiPath, from, to, dc);

		}

		// fixes or floating or ... ql 에서 지원되는 coupon들 여기서 계산
		Real legVanillaAmt(const Leg& leg, const ScenarioPath& scenPath, Date excerciseDate, Date next_excerciseDate, Size& couponCalculatedCount) const
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
				// paymentDate in (excerciseDate, next_excerciseDate]
				if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
					continue;

				couponCalculatedCount += 1;
				DayCounter dc = fltCpn->dayCounter();

				Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
				Real payment_t = dc.yearFraction(evaluationDate, paymentDate);
				Real reset_t = dc.yearFraction(evaluationDate, fltCpn->fixingDate());

				Real df = this->discount(scenPath, excerciseDate, paymentDate, dc);
				// Real rate = this->rate(multiPath, ex_t);

				if (reset_t < 0)
				{
					amount += fltCpn->amount() * df;
				}
				else
				{


					//Real coupon_discounted = coupon * this->discount(multiPath, excerciseDate, paymentDate, dc); // 할인함
					Real coupon_discounted = coupon * df;
					amount += coupon_discounted;
				}
			}

			return amount;
		}

		// excerciseDate 시점의 continuation value
		// 최종적으로 이거는 StructuredLeg class 만들어서 거기로 보내야함.

		// structured 중간에 vanilla?s
		Real legCpnMcAmt(const Leg& leg, const ScenarioPath& scenPath, Date excerciseDate, Date next_excerciseDate, Size& couponCalculatedCount) const
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

			for (Size j = startLegCount; j < legNum; j++)
			{
				boost::shared_ptr<CouponMC> cpn = boost::dynamic_pointer_cast<CouponMC>(leg[j]);

				Date paymentDate = cpn->date();

				// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
				if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
					continue;

				couponCalculatedCount += 1;
				DayCounter dc = cpn->dayCounter();

				Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
				Real payment_t = dc.yearFraction(evaluationDate, paymentDate);

				// Real rate = this->rate(scenPath, ex_t);
				Real df = this->discount(scenPath, excerciseDate, paymentDate, dc);

				Real notional = cpn->notional();
				Real accrualPeriod = cpn->accrual();
				// 여기서 cms range를 계산한다고 쳐.

				Real coupon = notional * accrualPeriod * cpn->calculate_path(scenPath);
				Real coupon_discounted = coupon * this->discount(scenPath, excerciseDate, paymentDate, dc); // 할인함

				amount += coupon_discounted;
			}

			return amount;

		}


		//Real _cpnAmt_mc(const boost::shared_ptr<CouponMC>& cpn, const ScenarioPath& scenPath) const
		//{
		//	// Real notional = cpn->notional();
		//	// Real accrualPeriod = cpn->accrual();
		//	// 여기서 cms range를 계산한다고 쳐.

		//	// Real coupon = notional * accrualPeriod * cpn->calculate_path(scenPath);

		//	return cpn->calculate_path(scenPath);
		//}

		//Real _cpnAmt_flt(const boost::shared_ptr<FloatingRateCouponMC>& cpn, const ScenarioPath& scenPath) const
		//{
		//	Real amount = 0.0;
		//	DayCounter dc = cpn->dayCounter();

		//	// Real ex_t = dc.yearFraction(evaluationDate, excerciseDate);
		//	// Real payment_t = dc.yearFraction(evaluationDate, paymentDate);
		//	Real reset_t = dc.yearFraction(evaluationDate, cpn->fixingDate());

		//	// Real rate = this->rate(multiPath, ex_t);

		//	if (reset_t < 0)
		//	{
		//		amount = cpn->amount();
		//	}
		//	else
		//	{
		//		Real notional = cpn->notional();
		//		Real accrualPeriod = cpn->accrual();
		//		// 여기서 cms range를 계산한다고 쳐.

		//		Real amount = notional * accrualPeriod * cpn->calculate_path(scenPath);

		//	}

		//	return amount;

		//}
		//Real cpnAmt(const boost::shared_ptr<CashFlow>& cf, const ScenarioPath& scenPath) const
		//{
		//	Real v = 0.0;

		//	auto cpnMc = boost::dynamic_pointer_cast<CouponMC>(cf);

		//	cpnMc->calculate_path();

		//	if (boost::dynamic_pointer_cast<CouponMC>(cf))
		//	{
		//		v = _cpnAmt_mc(boost::dynamic_pointer_cast<CouponMC>(cf), scenPath);
		//	}
		//	else if (boost::dynamic_pointer_cast<FloatingRateCouponMC>(cf))
		//	{
		//		v = _cpnAmt_flt(boost::dynamic_pointer_cast<FloatingRateCouponMC>(cf), scenPath);

		//	}
		//	else if (boost::dynamic_pointer_cast<FixedRateCoupon>(cf))
		//	{

		//	}
		//	else
		//	{
		//		QL_FAIL("unknown cashflow");
		//	}

		//	return v;

		//}

		Real legAmt(const Leg& leg, const ScenarioPath& scenPath, Date excerciseDate, Date next_excerciseDate, Size& couponCalculatedCount, const DayCounter& dc) const
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

			for (Size j = startLegCount; j < legNum; j++)
			{
				auto cpn = boost::dynamic_pointer_cast<CouponMC>(leg[j]);

				Date paymentDate = cpn->date();

				// excerciseDate 와 next_excerciseDate 사이에 있는거만 넣음
				if (paymentDate <= excerciseDate || next_excerciseDate < paymentDate)
					continue;

				Real discount = this->discount(scenPath, excerciseDate, paymentDate, dc); // 할인함

				amount += cpn->calculate_path(scenPath) * discount;
			}

			return amount;

		}

		bool in_array(const std::string& value, const std::vector<std::string>& array) const
		{
			return std::find(array.begin(), array.end(), value) != array.end();
		}

		void calculate() const
		{
			const Leg& payLeg = this->arguments_.payLeg;
			const Leg& recLeg = this->arguments_.recLeg;

			// index check : coupon에서 사용된 index 가 전부 있는가?
#pragma region Check Index 
			std::vector<std::string> index_names;
			for (Size i = 0; i < this->indexes_.size(); i++)
				index_names.push_back(this->indexes_[i]->model()->name());

			for (Size i = 0; i < payLeg.size(); i++)
			{
				// check only case of CouponMC
				const boost::shared_ptr<CouponMC>& cpnMC = boost::dynamic_pointer_cast<CouponMC>(payLeg[i]);

				if (cpnMC)
				{
					std::vector<std::string> names = cpnMC->index_names();
					for (Size j = 0; j < names.size(); j++)
						QL_REQUIRE(in_array(names[j], index_names), names[j] << " does not exist in scenario");
				}
			}
#pragma endregion

			Real side = 1.0;

			Date evaluationDate = Settings::instance().evaluationDate();

			//std::vector<boost::function1<Real, Real>> v =
			std::vector<boost::function1<Real, Array>> v =
				LsmBasisSystem::multiPathBasisSystem(2, 2, LsmBasisSystem::Hermite);

			const Size n = this->samples_;
			std::vector<Real> prices(n), exercise(n);

			// precalculate!

			DayCounter dc = Actual365Fixed();
			DayCounter dc_360 = Thirty360();

			const TimeGrid& timeGrid = this->timeGrid_;

			// std::vector<Time> excerciseTimes;

			std::vector<Date> excerciseDates;
			for (Size i = 0; i < this->arguments_.excerciseDates.size(); i++)
			{
				if (evaluationDate < this->arguments_.excerciseDates[i])
					excerciseDates.push_back(this->arguments_.excerciseDates[i]);
			}

			Date maturityDate = payLeg.back()->date() >= recLeg.back()->date() ? payLeg.back()->date() : recLeg.back()->date();
			if (excerciseDates.back() != maturityDate)
				excerciseDates.push_back(maturityDate);

			const Size len = excerciseDates.size();
			coeff_.reset(new Array[len - 1]);

			std::vector<std::vector<Real>> precalculated_fixed_amts;
			std::vector<std::vector<Real>> precalculated_float_amts;
			std::vector<std::vector<Real>> precalculated_discounts;

			//std::vector<MultiPath> multiPaths;
			std::vector<ScenarioPath> scenPaths;
			std::vector<Real> floatingAmts;
			std::vector<Real> structuredAmts;
			std::vector<Real> dfs;
			//std::vector<Real> x;
			std::vector<Array> x;
			std::vector<Real> y;

			Size couponCalculatedCount = 0;
			Real floatLegNpv = 0.0;
			Real structuredLegNpv = 0.0;

			std::vector<Array> factors_arr;

			for (Size i = 0; i < this->indexes_.size(); i++)
				factors_arr.push_back(Array(this->indexes_[i]->model()->factors(), 0.0));

			// 여기서 Scenario를 gen을 다해.
			// 그리고 그 file을 이용하는 방법...? -> 무조건 file...? 이게 나을 거같은...
			// 그러면 두번 일하는 거도 없고, 그 시나리오를 통해서 계산하는 거도 자연스러움.
			// 파일 생성을 안한다면, 메모리에 넣어야하는데, 그거는...? 그냥 따로 짜야함.
			// 우선 먼저 file을 이용하는 거로 함.
			// 아니면 그냥 메모리에 죄다 박어.

			std::vector<boost::shared_ptr<ProcessValue>> models;
			for (Size i = 0; i < this->indexes_.size(); i++)
				models.push_back(this->indexes_[i]->model());

			const boost::shared_ptr<QuantLib::IRsgWrapper>& rsgWrapper = this->rsgWrapper_;

			bool autoGenerate = true;
			const std::string& filename = "dummy";
			bool processMomentMatch = false;

			// discount 넣어야하네... 
			// calcs
			// python.......?!

			// 두번 계산되는 데...? -> 외부에서 들어오는 경우에...? 

			std::vector<boost::shared_ptr<ProcessValue>> calcs;

			for (Size i = 0; i < indexes_.size(); i++)
			{
				calcs.push_back(indexes_[i]->get_additionalCalc());

			}

			boost::shared_ptr<DiscountFactorCalc> pay_discountFactor(new DiscountFactorCalc("pay_discount", this->pay_discount_model_));
			calcs.push_back(pay_discountFactor);
			boost::shared_ptr<DiscountFactorCalc> rec_discountFactor(new DiscountFactorCalc("rec_discount", this->rec_discount_model_));
			calcs.push_back(rec_discountFactor);

			EvolverFactory::scenario_generator2(models, calcs, this->corr_, this->timeGrid_, rsgWrapper, filename, processMomentMatch);
			ScenarioResultReader scenarioResult(filename);

			for (Size j = 0; j < n; ++j)
			{
				ScenarioPath scenPath(scenarioResult.multiPath(j), timeGrid_);
				scenPaths.push_back(scenPath);

				prices[j] = 0.0;

				Real float_amt = this->legFltAmt(floatingLeg, scenPath, evaluationDate, maturityDate, couponCalculatedCount);
				Real structured_amt = this->legStructuredAmt(structuredLeg, scenPath, evaluationDate, maturityDate, couponCalculatedCount);

				floatLegNpv += float_amt;
				structuredLegNpv += structured_amt;

				std::cout << j << " -> " << (floatLegNpv - structuredLegNpv) / (j + 1) << " = " << floatLegNpv / (j + 1) << " " << structuredLegNpv / (j + 1) << std::endl;
			}


			floatLegNpv = floatLegNpv / n;
			structuredLegNpv = structuredLegNpv / n;

			std::cout << floatLegNpv - structuredLegNpv << " " << floatLegNpv << " " << structuredLegNpv << std::endl;
			std::cout << "option pricing --------------------------------------------------" << std::endl;

			//for (int i = 2; i >= 0; --i) {
			for (int i = len - 2; i >= 0; --i) {
				x.clear();
				y.clear();
				dfs.clear();
				floatingAmts.clear();
				structuredAmts.clear();

				Time t = dc.yearFraction(evaluationDate, excerciseDates[i]);
				Size pos = this->timeGrid_.closestIndex_Date(excerciseDates[i]);
				//this->hw_model_ = boost::shared_ptr<HullWhite>(
				//	new HullWhite(this->process_->fitting(), this->process_->a(t), this->process_->sigma(t)));

				//roll back step
				for (Size j = 0; j < n; ++j)
				{
					const MultiPath& multiPath = scenPaths[j];

					// Real rate = this->rate(multiPath, t);
					Array factors(2);
					factors[0] = multiPath[0][pos];
					factors[1] = multiPath[1][pos];

					Real df = this->discount(multiPath, excerciseDates[i], excerciseDates[i + 1], dc);

					Real float_amt = this->legFltAmt(floatingLeg, multiPath, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);
					Real structured_amt = this->legStructuredAmt(structuredLeg, scenPath, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);

					dfs.push_back(df);
					floatingAmts.push_back(float_amt);
					structuredAmts.push_back(structured_amt);

					Real price = prices[j] * df + float_amt - structured_amt;
					Array states(2);

					states[0] = this->indexes_[0]->calc_index(t, factors, this->model_);
					states[1] = this->indexes_[1]->calc_index(t, factors, this->model_);

					x.push_back(states);
					y.push_back(price);
				}

				coeff_[i] = GeneralLinearLeastSquares(x, y, v).coefficients();

				for (Size j = 0; j < n; ++j) {
					const MultiPath& multiPath = multiPaths[j];
					//Real rate = this->rate(multiPath, t);
					Array factors(2);
					factors[0] = multiPath[0][pos];
					factors[1] = multiPath[1][pos];

					//Real float_amt = this->legFltAmt(floatingLeg, multiPath, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);
					//Real structured_amt = this->legStructuredAmt(structuredLeg, multiPath, excerciseDates[i], excerciseDates[i + 1], couponCalculatedCount);

					//Real df = this->discount(multiPath, excerciseDates[i], excerciseDates[i + 1], dc);

					Real float_amt = floatingAmts[j];
					Real structured_amt = structuredAmts[j];
					Real df = dfs[j];

					Real price = prices[j] * df + float_amt - structured_amt;
					prices[j] = price;
					exercise[j] = 0.0;

					Array states(2);

					states[0] = this->indexes_[0]->calc_index(t, factors, this->model_);
					states[1] = this->indexes_[1]->calc_index(t, factors, this->model_);

					Real continuationValue = 0.0;

					for (Size l = 0; l < v.size(); ++l) {
						continuationValue += coeff_[i][l] * v[l](states);
					}

					if (continuationValue < exercise[j]) {
						prices[j] = exercise[j];
					}
				}

				std::cout << i << " / " << len - 2 << std::endl;
			}

			// 평가일까지 한번 더 땡겨야함.
			for (Size j = 0; j < n; ++j) {

				const MultiPath& multiPath = multiPaths[j]; // 어차피 첫번째꺼 사용할 거임. 같은값나오므로.
				Real float_amt = this->legFltAmt(floatingLeg, multiPath, evaluationDate, excerciseDates[0], couponCalculatedCount);
				Real structured_amt = this->legStructuredAmt(structuredLeg, multiPath, evaluationDate, excerciseDates[0], couponCalculatedCount);

				Real df = this->discount(multiPath, evaluationDate, excerciseDates[0], dc);
				Real price = prices[j] * df + float_amt - structured_amt;
				prices[j] = price;
			}

			Real npv = 0.0;

			for (Size i = 0; i < n; i++)
			{
				npv += prices[i] / n;
			}

			this->results_.value = npv;

			std::cout << npv << std::endl;
		}

	private:
		Size samples_;
		Size polynomOrder_;
		// LsmBasisSystem::PolynomType polynomType_;
		mutable boost::scoped_array<Array> coeff_;

		Parameter phi_;

		Size timeSteps_;
		Size seed_;
		Size nCalibrationSamples_;
		TimeGrid timeGrid_;

		std::vector<boost::shared_ptr<IndexMC>> indexes_;
		QuantLib::Matrix corr_;
		boost::shared_ptr<ProcessModel> pay_discount_model_;
		boost::shared_ptr<ProcessModel> rec_discount_model_;


	};

}


namespace StructuredSwapDev
{
	Size simulNum()
	{
		return 10000;
	}

	boost::shared_ptr<SwapIndexMC> swap_index(const std::string& familyName, Period tenor,
		const boost::shared_ptr<IborIndex>& iborIndex,
		const boost::shared_ptr<AffineModel>& model,
		const TimeGrid& timeGrid)
	{
		Natural settlementDays = 1;
		Currency currency = KRWCurrency();
		Period fixedLegTenor(3, TimeUnit::Months);

		auto swapindex = boost::shared_ptr<SwapIndex>(new
			SwapIndex(familyName,
				tenor,
				settlementDays,
				currency,
				iborIndex->fixingCalendar(),
				fixedLegTenor,
				iborIndex->businessDayConvention(),
				iborIndex->dayCounter(),
				iborIndex));

		return boost::shared_ptr<SwapIndexMC>(new SwapIndexMC(swapindex, model, timeGrid));
	}

	boost::shared_ptr<YieldTermStructure> irs_curve()
	{
		Date evaluationDate = Settings::instance().evaluationDate();

		// 2021-12-31
		std::vector<MarketCurveRate> marketCurveRates
		{
			MarketCurveRate("1D", 0.01436, MarketCurveRate::Type::Cash),
			MarketCurveRate("3M", 0.0129, MarketCurveRate::Type::Cash),
			MarketCurveRate("6M", 0.01405, MarketCurveRate::Type::Swap),
			MarketCurveRate("9M", 0.01485, MarketCurveRate::Type::Swap),
			MarketCurveRate("1Y", 0.015825, MarketCurveRate::Type::Swap),
			MarketCurveRate("18M", 0.01705, MarketCurveRate::Type::Swap),
			MarketCurveRate("2Y", 0.0177, MarketCurveRate::Type::Swap),
			MarketCurveRate("3Y", 0.0182, MarketCurveRate::Type::Swap),
			MarketCurveRate("4Y", 0.018375, MarketCurveRate::Type::Swap),
			MarketCurveRate("5Y", 0.018375, MarketCurveRate::Type::Swap),
			MarketCurveRate("6Y", 0.0183, MarketCurveRate::Type::Swap),
			MarketCurveRate("7Y", 0.01835, MarketCurveRate::Type::Swap),
			MarketCurveRate("8Y", 0.018375, MarketCurveRate::Type::Swap),
			MarketCurveRate("9Y", 0.018475, MarketCurveRate::Type::Swap),
			MarketCurveRate("10Y", 0.01855, MarketCurveRate::Type::Swap),
			MarketCurveRate("12Y", 0.018475, MarketCurveRate::Type::Swap),
			MarketCurveRate("15Y", 0.0174, MarketCurveRate::Type::Swap),
			MarketCurveRate("20Y", 0.0162, MarketCurveRate::Type::Swap),
			MarketCurveRate("25Y", 0.016275, MarketCurveRate::Type::Swap),
			MarketCurveRate("30Y", 0.01625, MarketCurveRate::Type::Swap)
		};

		return YieldCurveExt::bootstrapping_ccp(evaluationDate, marketCurveRates,
			Interpolator1D::Linear,
			Extrapolator1D::FlatForward, "IRSKRW_KRCCP");
	}

	boost::shared_ptr<G2ExtProcess> g2ext_process(Real sigma1, Real sigma2)
	{
		std::vector<Real> alpha_times{ 100.0 };

		PiecewiseConstantParameter alpha1_para(alpha_times);
		alpha1_para.setParam(0, 0.05);
		PiecewiseConstantParameter alpha2_para(alpha_times);
		alpha2_para.setParam(0, 0.5);

		std::vector<Real> sigma_times{ 100.0 };
		PiecewiseConstantParameter sigma1_para(sigma_times);
		sigma1_para.setParam(0, sigma1);
		PiecewiseConstantParameter sigma2_para(sigma_times);
		sigma2_para.setParam(0, sigma2);

		Real corr = -0.9;

		//boost::shared_ptr<HullWhite1FactorProcess> process(
		//	new HullWhite1FactorProcess(Handle<YieldTermStructure>(irs_curve()), alpha_para, sigma_para));
		boost::shared_ptr<G2ExtProcess> process(
			new G2ExtProcess(Handle<YieldTermStructure>(irs_curve()), alpha1_para, sigma1_para, alpha2_para, sigma2_para, corr));

		return process;

	}

	// dual spread
	void test_24537()
	{
		Date effectiveDate = DateParser::parseISO("2017-01-17");

		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 20000000000;

		Leg structuredLeg;

		Real upperTrigger = 0.055;
		Real inCoupon = 0.0356;

		std::vector<Date> structured_payment_dates
		{
			DateParser::parseISO("2018-01-17"),
			DateParser::parseISO("2019-01-17"),
			DateParser::parseISO("2020-01-17"),
			DateParser::parseISO("2021-01-18"),
			DateParser::parseISO("2022-01-17"), // 2.955 -> 2.877342
			DateParser::parseISO("2023-01-17"),
			DateParser::parseISO("2024-01-17"),
			DateParser::parseISO("2025-01-17"),
			DateParser::parseISO("2026-01-19"),
			DateParser::parseISO("2027-01-18"),
			DateParser::parseISO("2028-01-17"),
			DateParser::parseISO("2029-01-17"),
			DateParser::parseISO("2030-01-17"),
			DateParser::parseISO("2031-01-17"),
			DateParser::parseISO("2032-01-19"),

		};

		DayCounter dayCounter = Actual365Fixed();
		Date maturityDate = structured_payment_dates.back();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		//TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		auto process = g2ext_process(0.007287314, 0.002576387);

		const boost::shared_ptr<AffineModel>& g2 =
			boost::dynamic_pointer_cast<AffineModel>(boost::shared_ptr<G2>(
				new G2(process->termStructure(),
					process->a(0.0), process->sigma(0.0),
					process->b(0.0), process->eta(0.0),
					process->rho())));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, g2, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS5Y", Period(5, TimeUnit::Years), index, g2, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, g2, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 4) ? 0.02877342 : 0.0;  // pre accruedAmount

			boost::shared_ptr<DualSpreadRangeAccrualCouponMC> cpn(
				new DualSpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					indexMC3,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					upperTrigger,
					inCoupon, 0.0, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Leg floatingLeg;

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.005;

		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-04-17"), notional, DateParser::parseISO("2017-01-17"), DateParser::parseISO("2017-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-07-17"), notional, DateParser::parseISO("2017-04-17"), DateParser::parseISO("2017-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2017-10-17"), notional, DateParser::parseISO("2017-07-17"), DateParser::parseISO("2017-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-01-17"), notional, DateParser::parseISO("2017-10-17"), DateParser::parseISO("2018-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-04-17"), notional, DateParser::parseISO("2018-01-17"), DateParser::parseISO("2018-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-07-17"), notional, DateParser::parseISO("2018-04-17"), DateParser::parseISO("2018-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2018-10-17"), notional, DateParser::parseISO("2018-07-17"), DateParser::parseISO("2018-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-01-17"), notional, DateParser::parseISO("2018-10-17"), DateParser::parseISO("2019-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-04-17"), notional, DateParser::parseISO("2019-01-17"), DateParser::parseISO("2019-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-07-17"), notional, DateParser::parseISO("2019-04-17"), DateParser::parseISO("2019-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2019-10-17"), notional, DateParser::parseISO("2019-07-17"), DateParser::parseISO("2019-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-01-17"), notional, DateParser::parseISO("2019-10-17"), DateParser::parseISO("2020-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-04-17"), notional, DateParser::parseISO("2020-01-17"), DateParser::parseISO("2020-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-07-17"), notional, DateParser::parseISO("2020-04-17"), DateParser::parseISO("2020-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2020-10-19"), notional, DateParser::parseISO("2020-07-17"), DateParser::parseISO("2020-10-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-01-18"), notional, DateParser::parseISO("2020-10-19"), DateParser::parseISO("2021-01-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-04-19"), notional, DateParser::parseISO("2021-01-18"), DateParser::parseISO("2021-04-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-07-19"), notional, DateParser::parseISO("2021-04-19"), DateParser::parseISO("2021-07-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2021-10-18"), notional, DateParser::parseISO("2021-07-19"), DateParser::parseISO("2021-10-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-01-17"), notional, DateParser::parseISO("2021-10-18"), DateParser::parseISO("2022-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-04-18"), notional, DateParser::parseISO("2022-01-17"), DateParser::parseISO("2022-04-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-07-18"), notional, DateParser::parseISO("2022-04-18"), DateParser::parseISO("2022-07-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2022-10-17"), notional, DateParser::parseISO("2022-07-18"), DateParser::parseISO("2022-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-01-17"), notional, DateParser::parseISO("2022-10-17"), DateParser::parseISO("2023-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-04-17"), notional, DateParser::parseISO("2023-01-17"), DateParser::parseISO("2023-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-07-17"), notional, DateParser::parseISO("2023-04-17"), DateParser::parseISO("2023-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2023-10-17"), notional, DateParser::parseISO("2023-07-17"), DateParser::parseISO("2023-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-01-17"), notional, DateParser::parseISO("2023-10-17"), DateParser::parseISO("2024-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-04-17"), notional, DateParser::parseISO("2024-01-17"), DateParser::parseISO("2024-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-07-17"), notional, DateParser::parseISO("2024-04-17"), DateParser::parseISO("2024-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2024-10-17"), notional, DateParser::parseISO("2024-07-17"), DateParser::parseISO("2024-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-01-17"), notional, DateParser::parseISO("2024-10-17"), DateParser::parseISO("2025-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-04-17"), notional, DateParser::parseISO("2025-01-17"), DateParser::parseISO("2025-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-07-17"), notional, DateParser::parseISO("2025-04-17"), DateParser::parseISO("2025-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2025-10-17"), notional, DateParser::parseISO("2025-07-17"), DateParser::parseISO("2025-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-01-19"), notional, DateParser::parseISO("2025-10-17"), DateParser::parseISO("2026-01-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-04-17"), notional, DateParser::parseISO("2026-01-19"), DateParser::parseISO("2026-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-07-17"), notional, DateParser::parseISO("2026-04-17"), DateParser::parseISO("2026-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2026-10-19"), notional, DateParser::parseISO("2026-07-17"), DateParser::parseISO("2026-10-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-01-18"), notional, DateParser::parseISO("2026-10-19"), DateParser::parseISO("2027-01-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-04-19"), notional, DateParser::parseISO("2027-01-18"), DateParser::parseISO("2027-04-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-07-19"), notional, DateParser::parseISO("2027-04-19"), DateParser::parseISO("2027-07-19"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2027-10-18"), notional, DateParser::parseISO("2027-07-19"), DateParser::parseISO("2027-10-18"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-01-17"), notional, DateParser::parseISO("2027-10-18"), DateParser::parseISO("2028-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-04-17"), notional, DateParser::parseISO("2028-01-17"), DateParser::parseISO("2028-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-07-17"), notional, DateParser::parseISO("2028-04-17"), DateParser::parseISO("2028-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2028-10-17"), notional, DateParser::parseISO("2028-07-17"), DateParser::parseISO("2028-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-01-17"), notional, DateParser::parseISO("2028-10-17"), DateParser::parseISO("2029-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-04-17"), notional, DateParser::parseISO("2029-01-17"), DateParser::parseISO("2029-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-07-17"), notional, DateParser::parseISO("2029-04-17"), DateParser::parseISO("2029-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2029-10-17"), notional, DateParser::parseISO("2029-07-17"), DateParser::parseISO("2029-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-01-17"), notional, DateParser::parseISO("2029-10-17"), DateParser::parseISO("2030-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-04-17"), notional, DateParser::parseISO("2030-01-17"), DateParser::parseISO("2030-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-07-17"), notional, DateParser::parseISO("2030-04-17"), DateParser::parseISO("2030-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2030-10-17"), notional, DateParser::parseISO("2030-07-17"), DateParser::parseISO("2030-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-01-17"), notional, DateParser::parseISO("2030-10-17"), DateParser::parseISO("2031-01-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-04-17"), notional, DateParser::parseISO("2031-01-17"), DateParser::parseISO("2031-04-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-07-17"), notional, DateParser::parseISO("2031-04-17"), DateParser::parseISO("2031-07-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2031-10-17"), notional, DateParser::parseISO("2031-07-17"), DateParser::parseISO("2031-10-17"), fixingDays, index, gearing, spread)));
		floatingLeg.push_back(boost::shared_ptr<IborCoupon>(new IborCoupon(DateParser::parseISO("2032-01-19"), notional, DateParser::parseISO("2031-10-17"), DateParser::parseISO("2032-01-19"), fixingDays, index, gearing, spread)));



		std::vector<Date> optioin_exdates
		{
			DateParser::parseISO("2019-01-17"),
			DateParser::parseISO("2020-01-17"),
			DateParser::parseISO("2021-01-18"),
			DateParser::parseISO("2022-01-17"),
			DateParser::parseISO("2023-01-17"),
			DateParser::parseISO("2024-01-17"),
			DateParser::parseISO("2025-01-17"),
			DateParser::parseISO("2026-01-19"),
			DateParser::parseISO("2027-01-18"),
			DateParser::parseISO("2028-01-17"),
			DateParser::parseISO("2029-01-17"),
			DateParser::parseISO("2030-01-17"),
			DateParser::parseISO("2031-01-17"),
		};

		index->addFixing(DateParser::parseISO("2021-10-15"), 0.0107);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007287314, 0.002576387),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_31664()
	{
		Date effectiveDate = DateParser::parseISO("2016-11-17");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2031-11-17");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 30000000000;

		Leg structuredLeg;

		Real upperTrigger = 0.06;
		Real inCoupon = 0.0311;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		//std::vector<Date> structured_payment_dates
		//{
		//	DateParser::parseISO("2018-01-17"),
		//	DateParser::parseISO("2019-01-17"),
		//	DateParser::parseISO("2020-01-17"),
		//	DateParser::parseISO("2021-01-18"),
		//	DateParser::parseISO("2022-01-17"), // 1.47
		//	DateParser::parseISO("2023-01-17"),
		//	DateParser::parseISO("2024-01-17"),
		//	DateParser::parseISO("2025-01-17"),
		//	DateParser::parseISO("2026-01-19"),
		//	DateParser::parseISO("2027-01-18"),
		//	DateParser::parseISO("2028-01-17"),
		//	DateParser::parseISO("2029-01-17"),
		//	DateParser::parseISO("2030-01-17"),
		//	DateParser::parseISO("2031-01-17"),
		//	DateParser::parseISO("2032-01-19"),

		//};

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS2Y", Period(2, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 20) ? 0.0147 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<DualSpreadRangeAccrualCouponMC> cpn(
				new DualSpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					indexMC3,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					upperTrigger,
					inCoupon, 0.0, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.0035;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 4; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-11-16"), 0.0115);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007397706801275404, 0.0025510499079574585),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_16651()
	{
		Date effectiveDate = DateParser::parseISO("2016-09-20");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2031-09-20");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 20000000000;

		Leg structuredLeg;

		Real upperTrigger = 0.055;
		Real inCoupon = 0.027;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(1, Years), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS5Y", Period(5, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 5) ? 0.00425914 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<DualSpreadRangeAccrualCouponMC> cpn(
				new DualSpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					indexMC3,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					upperTrigger,
					inCoupon, 0.0, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = -0.001;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 5; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-12-17"), 0.0127);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007277431759404901, 0.0024482344487102383),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_11728()
	{
		Date effectiveDate = DateParser::parseISO("2016-09-09");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2031-09-09");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 10000000000;

		Leg structuredLeg;

		Real upperTrigger = 0.06;
		Real inCoupon = 0.0305;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS5Y", Period(5, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 21) ? 0.00135339 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<DualSpreadRangeAccrualCouponMC> cpn(
				new DualSpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					indexMC3,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					upperTrigger,
					inCoupon, 0.0, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.0011;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 21; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-12-08"), 0.0127);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007277422768895635, 0.00231257708401405),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_01285()
	{
		Date effectiveDate = DateParser::parseISO("2016-09-05");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2031-09-05");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 5000000000;

		Leg structuredLeg;

		Real upperTrigger = 0.055;
		Real inCoupon = 0.0304;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS5Y", Period(5, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 21) ? 0.001484626 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<DualSpreadRangeAccrualCouponMC> cpn(
				new DualSpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					indexMC3,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					upperTrigger,
					inCoupon, 0.0, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.0008;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 21; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-12-03"), 0.0127);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007275954783768923, 0.0023018531246268062),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_13273()
	{
		Date effectiveDate = DateParser::parseISO("2017-02-17");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2032-02-17");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 50000000000;

		Leg structuredLeg;

		Real inCoupon = 0.0362;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(1, Years), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS2Y", Period(2, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 4) ? 0.03135 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<SpreadRangeAccrualCouponMC> cpn(
				new SpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					inCoupon, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.004;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 4; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-11-16"), 0.0115);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007392794838356776, 0.002722110685745046),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

	void test_16014()
	{
		Date effectiveDate = DateParser::parseISO("2016-06-17");
		// Date maturityDate = SouthKorea().advance(effectiveDate, Period(15, Years));
		Date maturityDate = DateParser::parseISO("2031-06-17");

		Calendar calendar = SouthKorea();
		Settings::instance().evaluationDate() = DateParser::parseISO("2021-12-31");
		Date evaluationDate = Settings::instance().evaluationDate();
		Real notional = 10000000000;

		Leg structuredLeg;

		Real inCoupon = 0.027;

		Schedule structured_schedule(effectiveDate, maturityDate, Period(1, Years), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		std::vector<Date> structured_payment_dates;

		for (Size i = 1; i < structured_schedule.size(); i++)
		{
			structured_payment_dates.push_back(
				structured_schedule.date(i));
		}

		DayCounter dayCounter = Actual365Fixed();

		//Size timeSteps = static_cast<Size>(Actual365Fixed().yearFraction(evaluationDate, maturityDate) * 365);
		//Time maturityTime = dayCounter.yearFraction(evaluationDate, structured_payment_dates.back());
		//TimeGrid timeGrid(maturityTime, timeSteps);

		// TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 365, "CUSTOM", 1, 1);
		TimeGrid timeGrid = TestEngine::build_timeGrid(evaluationDate, 11, 52, "CUSTOM", 1, 1);

		for (Size i = 0; i < timeGrid.size() - 1; i++)
		{
			if (timeGrid.date_at(i) == timeGrid.date_at(i + 1))
				std::cout << timeGrid.date_at(i) << std::endl;
		}

		Handle<YieldTermStructure> yts_h(irs_curve());

		boost::shared_ptr<IborIndex> index(
			new KrwCD(Period(3, TimeUnit::Months), yts_h));

		const boost::shared_ptr<SwapIndexMC>& indexMC1 = swap_index("CMS10Y", Period(10, TimeUnit::Years), index, timeGrid);
		const boost::shared_ptr<SwapIndexMC>& indexMC2 = swap_index("CMS5Y", Period(5, TimeUnit::Years), index, timeGrid);

		boost::shared_ptr<IborIndexMC> indexMC3(new IborIndexMC(index, timeGrid));

		for (Size i = 0; i < structured_payment_dates.size(); i++)
		{
			Date accrualStartDate = i == 0 ? effectiveDate : structured_payment_dates[i - 1];
			Date accrualEndDate = structured_payment_dates[i];

			Real accruedAmount = (i == 5) ? 0.010333 : 0.0;  // pre accruedAmount
			// Real accruedAmount = (i == 20) ? 0.0 : 0.0;  // for clean -> 이거 이렇게 돌리면 안됨. 나중에 옵션계산할때 dirty 기준으로 해야함. 
			// accured를 나중에 계산해야함. // 

			boost::shared_ptr<SpreadRangeAccrualCouponMC> cpn(
				new SpreadRangeAccrualCouponMC(
					structured_payment_dates[i],
					notional,
					indexMC1, //
					indexMC2,
					accrualStartDate,
					accrualEndDate,
					dayCounter,
					timeGrid,
					inCoupon, 0.0, accruedAmount));

			structuredLeg.push_back(cpn);
		}

		Natural fixingDays = 1;
		Real gearing = 1.0;
		Real spread = 0.001;


		Schedule floating_schedule(effectiveDate, maturityDate, Period(3, Months), calendar, ModifiedFollowing, ModifiedFollowing, DateGeneration::Backward, false);
		auto floatingLeg = IborLeg(floating_schedule, index)
			.withNotionals(notional)
			.withPaymentDayCounter(dayCounter)
			.withSpreads(spread);


		std::vector<Date> optioin_exdates;
		for (Size i = 5; i < structured_schedule.size() - 1; i++)
		{
			optioin_exdates.push_back(
				structured_schedule.date(i));
		}

		index->addFixing(DateParser::parseISO("2021-12-16"), 0.0127);

		boost::shared_ptr<LegExerciseOption> option(
			new LegExerciseOption(optioin_exdates, structured_payment_dates));

		boost::shared_ptr<Swap> swap(
			new TestEngine::StructuredSwap(structuredLeg, floatingLeg, option));

		Size samples = simulNum();

		std::vector<boost::shared_ptr<SwapIndexMC>> indexes{ indexMC1, indexMC2 };
		boost::shared_ptr<MCStructuredSwapType1Engine> engine(
			new MCStructuredSwapType1Engine(
				g2ext_process(0.007272434441901232, 0.002396455905339247),
				timeGrid,
				samples,
				1,
				2048,
				indexes));

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

}