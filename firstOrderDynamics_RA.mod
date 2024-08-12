// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials)
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

// Define economic parameters
parameters bbeta ssigma aalpha ddelta aggEmployment
	mmu ttau rrhoTFP ssigmaTFP;
bbeta = .96;										% discount factor (annual calibration)
ssigma = 1;											% coefficient of relative risk aversion
aalpha = .36;										% capital share
ddelta = .1;										% depreciation rate (annual calibration)
aggEmployment = .93;
mmu = .15;
ttau = mmu * (1 - aggEmployment) / aggEmployment;
rrhoTFP = .859;										
ssigmaTFP = .014;
//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var r w;

//----------------------------------------------------------------
// asset captial
//----------------------------------------------------------------

var assets aggregateCapital;

//----------------------------------------------------------------
// Aggregate TFP
//----------------------------------------------------------------

var aggregateTFP;

//----------------------------------------------------------------
// Auxiliary variables we're interested in
//----------------------------------------------------------------

var logAggregateOutput logAggregateInvestment logAggregateConsumption logWage;

//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

varexo aggregateTFPShock;

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

//----------------------------------------------------------------
// Factor prices (# equations = 2)
//----------------------------------------------------------------

aggregateCapital = assets(-1);
r = exp(aggregateTFP) * aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

//----------------------------------------------------------------
// Law of motion for aggregate TFP (# equations = 1)
//----------------------------------------------------------------

aggregateTFP = rrhoTFP * aggregateTFP(-1) + ssigmaTFP * aggregateTFPShock;

//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 4)
//----------------------------------------------------------------

// Output
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));

// Investment
logAggregateInvestment = log(assets - (1 - ddelta) * aggregateCapital);

// Consumption
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));

// Eulear
(1+r)*bbeta * exp(ssigma*(logAggregateConsumption-logAggregateConsumption(+1)))=1;

// Wage
logWage = log(w);

end;

//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
end;

steady;

// Simulate
stoch_simul(order=1,hp_filter=100,irf=40) aggregateTFP logAggregateOutput 
	logAggregateConsumption logAggregateInvestment logWage r;
