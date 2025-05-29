// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials)
//
// Thomas Winberry, July 26th, 2016

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

// Define economic parameters
parameters bbeta ssigma aalpha ddelta efficiency
	rrhoTFP ssigmaTFP vvarphi cchi;
bbeta = .96;										% discount factor (annual calibration)
ssigma = 3;											% coefficient of relative risk aversion
aalpha = .36;										% capital share
ddelta = .1; % depreciation rate (annual calibration)
vvarphi = 1;
efficiency = .93;
rrhoTFP = .859;										
ssigmaTFP = .014;
ASS = 1;
lss = 1;
rss = 1/bbeta-1;
K_L = (aalpha/(rss+ddelta))^(1/(1-aalpha));
KSS = K_L * efficiency * lss;
ISS = ddelta * KSS;
YSS = Y_L * efficiency * lss;
Y_L = ASS * K_L^aalpha;
CSS = YSS- ISS;
wss = (1-aalpha) * K_L^aalpha;
cchi = lss^vvarphi * CSS^ssigma / wss;

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

//----------------------------------------------------------------
// Prices
//----------------------------------------------------------------

var r w;

//----------------------------------------------------------------
// labor asset captial
//----------------------------------------------------------------

var labor assets aggregateCapital;

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
r = exp(aggregateTFP) * aalpha * (aggregateCapital ^ (aalpha - 1)) * ((efficiency*labor) ^ (1 - aalpha)) - ddelta;
w = exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (1 - aalpha) * ((efficiency*labor) ^ (-aalpha));

//----------------------------------------------------------------
// Law of motion for aggregate TFP (# equations = 1)
//----------------------------------------------------------------

aggregateTFP = rrhoTFP * aggregateTFP(-1) + ssigmaTFP * aggregateTFPShock;

//----------------------------------------------------------------
// Auxiliary variables of interest (# equations = 4)
//----------------------------------------------------------------

// Output
logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * ((efficiency*labor) ^ (1 - aalpha)));

// Investment
logAggregateInvestment = log(assets - (1 - ddelta) * aggregateCapital);

// Consumption
logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));

// Eulear
(1+r)* bbeta * exp(ssigma*(logAggregateConsumption-logAggregateConsumption(+1)))=1;

// labor supply
cchi*w=labor^vvarphi*exp(ssigma*logAggregateConsumption);
// Wage
logWage = log(w);

end;


//----------------------------------------------------------------
// 4. Computation
//----------------------------------------------------------------

initval;
    labor = lss;
    aggregateTFP = 0;
    r = rss;
    w = wss;
    logWage = log(w);
    aggregateCapital = KSS;
    assets = aggregateCapital;
    logAggregateOutput = log(YSS);
    logAggregateConsumption = log(CSS);
    logAggregateInvestment = log(ISS);
end;

steady;

// Specify shock process

shocks;
    var aggregateTFPShock = 1;
end;


// Simulate
// stoch_simul(order=1,hp_filter=100,irf=40);
stoch_simul(order=1,hp_filter=100,irf=40) aggregateTFP logAggregateOutput 
	logAggregateConsumption logAggregateInvestment labor r;
