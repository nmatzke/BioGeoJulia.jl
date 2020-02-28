



class Parameter {
	constructor(distNum, id, unicode = id, defaultVal = 1, minVal = null, maxVal = null, step = 0.01) {
		this.distNum = distNum;
		this.id = id;
		this.unicode = unicode;
		this.minVal = minVal;
		this.maxVal = maxVal;
		this.step = step;
		this.value = defaultVal;
		this.inclusiveLimit = [false, false];
	}
	
	getName() {
		return this.id;
	}

	loadFromDOM() {
		var span = $(`[distNum="` + this.distNum + `"][param="` + this.id + `"]`);
		if (span.length == 0 || span == null) return;
		var x = parseFloat(span.find("input").val());
		this.value = this.validate(x);
		span.find("input").val(this.value);
	}

	// Set a lower (0) or upper (1) limit to inclusive instead of exclusive (default)
	setInclusive(limit, invlusive = true) {
		this.inclusiveLimit[limit] = invlusive;
	}	

	set(x) {
		this.value = this.validate(x);
		console.log("Setting", this.id, "to", this.value);
	}

	get() {
		return this.value;
	}

	getHTML() {
		return `<span class="distr" distNum="` + this.distNum + `" param="` + this.id + `">
					` + this.unicode + ` = <input type="number" step="` + this.step + `" min="` + this.minVal + `" max = "` + this.maxVal + `" class="numberinput cranberry" onclick="$(this).select()" style="width:5em;" value="` + this.value + `" onchange="setParams();"> 
			 </span>`;
	}

	validate(x) {
		if (x == null) return 1; 
		if (this.minVal != null) if (x < this.minVal || (!this.inclusiveLimit[0] && x == this.minVal)) return this.minVal + 0.5 * this.step;
		if (this.maxVal != null) if (x > this.maxVal || (!this.inclusiveLimit[1] && x == this.maxVal)) return this.maxVal - 0.5 * this.step;
		return x;
	}

	getUnicode() {
		return this.unicode;
	}

	getStep() {
		return this.step;
	}

	getMinVal() {
		return this.minVal;
	}

	getMaxVal() {
		return this.maxVal;
	}

	toJSON() {
		var JSON = {distNum: this.distNum, id : this.id, unicode: this.unicode, step: this.step, value: this.value};
		if (this.minVal != null) JSON.minVal = this.minVal;
		if (this.maxVal != null) JSON.maxVal = this.maxVal;
		return JSON;
	}

}



class Distribution {


	constructor(distNum, col = null) {
		
		this.distNum = distNum;
		this.col = col == null ? getDefaultColour(distNum) : col;
		this.parameters = [];
		this.desc = "";
		this.initialise();
		this.loadParams();
		
	}


	getName() {
		return this.constructor.name;
	}


	getParamHTML() {

		var html = `<div distr="` + this.getName() + `" distNum="` + this.distNum + `">`;
		for (var i = 0; i < this.parameters.length; i ++ ){
			html += this.parameters[i].getHTML();
		}
		return html + "</div>";
	}

	// Called before rendering plot. Loads parameters from DOM
	loadParams() {
		for (var i = 0; i < this.parameters.length; i ++ ){
			this.parameters[i].loadFromDOM();
		}
		this.prepareDensity();
	}

	// Precompute any terms before getDensity is called eg. normalisations
	prepareDensity() {

	}

	// Returns the probability density of x
	getDensity(x) {
		if (x <= 0 || x >= 1) return 0;
		return 1;
	}
	
	getDesc() {
		return this.desc;
	}

	setColour(col) {
		this.col = col;
	}

	getColour() {
		return this.col;
	}

	getXRange() {
		return [0,1];
	}

	// Return the maximum probability density over the x range  
	getYMax() {
		var ymax = 0;
		var xmin = this.getXRange()[0];
		var xmax = this.getXRange()[1];
		for (var i = 1; i < NUM_PIECES; i ++){
			var x = xmin + i*(xmax - xmin)/NUM_PIECES; 
			var y = this.getDensity(x);
			ymax = Math.max(ymax, y);
		}
		return ymax;
	}


	toJSON() {
		var JSON = {name: this.constructor.name, distNum: this.distNum, col: this.col};
		var paramJSON = [];
		for (var i = 0; i < this.parameters.length; i ++ ){
			paramJSON[i] = this.parameters[i].toJSON();
		}
		JSON.parameters = paramJSON;
		return JSON;
	}


	// Parse from a JSON in the format:
	// {p: [val1, val2, ...], col: col}
	// where val1, val2,.. are the parameter values in the same order specified in the class
	parseFromJSON(JSON) {

		// Set colour
		if (JSON.col != null) this.col = JSON.col;

		// Set the parameters
		for (var i = 0; i < this.parameters.length; i ++){
			if (i < JSON.p.length) this.parameters[i].set(parseFloat(JSON.p[i])); 
		}

	}


};



class Uniform extends Distribution {

	initialise() {
		this.lower = new Parameter(this.distNum, "lower", "lower", 0, null, null, 1);
		this.upper = new Parameter(this.distNum, "upper", "upper", 1, null, null, 1);
		this.parameters = [this.lower, this.upper];
		this.desc = "x is uniformly distributed between the lower and upper limits.";
	}


	getXRange() {
		return [this.lower.get(), this.upper.get()];
	}


	getDensity(x) {
		var lower = this.lower.get();
		var upper = this.upper.get();
		if (lower >= upper) return 0;
		if (x <= lower || x >= upper) return 0;	
		return 1 / (upper - lower);
	}

}



class Bactrian extends Distribution {

	initialise() {
		this.m = new Parameter(this.distNum, "m", "m", 0.9, 0, 1, 0.02);
		this.m.setInclusive(0);
		this.parameters = [this.m];
		this.desc = "A bimodal normal distribution with a mean of 0 and a standard deviation of 1-m^2.";
	}


	getXRange() {
		return [-2,2];
	}


	getDensity(x) {
		var m = this.m.get();
		var density = 1 / (2*Math.sqrt(2*Math.PI*(1-m*m)));
		density = density * (Math.exp( -(x + m)*(x + m) / (2 * (1-m*m)) ) + Math.exp( -(x - m)*(x - m) / (2 * (1-m*m)) ));
		return density;
	}

}



class Normal extends Distribution {

	initialise() {
		this.mu = new Parameter(this.distNum, "mu", "&mu;", 0, null, null, 0.5);
		this.sigma = new Parameter(this.distNum, "sigma", "&sigma;", 1, 0, null, 0.5);
		this.parameters = [this.mu, this.sigma];
		this.desc = "Normal distribution with mean mu and standard deviation sigma.";
	
	}


	getXRange() {
		var mu = this.mu.get();
		var sigma = this.sigma.get();
		return [mu - 2.5* sigma, mu + 2.5* sigma];
	}


	getDensity(x) {
		var mu = this.mu.get();
		var sigma = this.sigma.get();
		var density = 1 / (sigma * Math.sqrt(2*Math.PI)) * Math.exp(-1/2 * ((x-mu)/(sigma)) **2 );
		return density;
	}

}


class LogNormal extends Distribution {

	initialise() {
		this.mu = new Parameter(this.distNum, "mu", "&mu;", 0, null, null, 0.5);
		this.sigma = new Parameter(this.distNum, "sigma", "&sigma;", 1, 0, null, 0.25);
		this.parameters = [this.mu, this.sigma];
		this.desc = "Log-normal distribution where the mean and standard deviation are mu and sigma in log-space.";
	}


	getXRange() {
		var mu = this.mu.get();
		var sigma = this.sigma.get();
		return [0, Math.exp(mu) * 5];
	}


	getDensity(x) {
		if (x <= 0) return 0;
		var mu = this.mu.get();
		var sigma = this.sigma.get();
		var exponent =  ((Math.log(x)-mu)**2) /(2*sigma*sigma);
		var density = 1 / (x * sigma * Math.sqrt(2*Math.PI)) * Math.exp(-exponent);
		return density;
	}

}


class LogNormalMean extends Distribution {

	initialise() {
		this.mean = new Parameter(this.distNum, "mean", "mean", 1, 0, null, 0.5);
		this.sigma = new Parameter(this.distNum, "sigma", "&sigma;", 0.5, 0, null, 0.05);
		this.parameters = [this.mean, this.sigma];
		this.desc = "Log-normal distribution where the mean is 'mean' and the standard deviation in log-space is 'sigma'.";
	
	}


	getXRange() {
		var mean = this.mean.get();
		var sigma = this.sigma.get();
		var mu = Math.log(mean) - sigma*sigma/2;
		return [0, Math.exp(mu) * 3];
	}


	getDensity(x) {
		if (x <= 0) return 0;
		var mean = this.mean.get();
		var sigma = this.sigma.get();
		var mu = Math.log(mean) - sigma*sigma/2;
		var exponent =  ((Math.log(x)-mu)**2) /(2*sigma*sigma);
		var density = 1 / (x * sigma * Math.sqrt(2*Math.PI)) * Math.exp(-exponent);
		return density;
	}

}


class Beta extends Distribution {

	initialise() {
		this.alpha = new Parameter(this.distNum, "alpha", "&alpha;", 1, 0, null, 0.2);
		this.beta = new Parameter(this.distNum, "beta", "&beta;", 1, 0, null, 0.2);
		this.parameters = [this.alpha, this.beta];
		this.norm = 1;
		this.desc = "Beta distribution with shape alpha and scale beta.";
	
	}


	prepareDensity() {
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		this.norm = math.gamma(alpha) * math.gamma(beta) / math.gamma(alpha + beta);
	}

	getDensity(x) {
		if (x <= 0 || x >= 1) return 0;
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		var density = Math.exp((alpha-1)*Math.log(x) + (beta-1)*Math.log(1-x)) / this.norm;
		return density;
	}

}




class Exponential extends Distribution {

	initialise() {
		this.rate = new Parameter(this.distNum, "rate", "&lambda;", 1, 0, null, 0.2);
		this.parameters = [this.rate];
		this.desc = "Exponential distribution with a mean of 1/lambda.";
	}


	getXRange() {
		var rate = this.rate.get();
		return [0, Math.log(2) / rate * 5];
	}

	getDensity(x) {
		if (x <= 0) return 0;
		var rate = this.rate.get();
		var density = rate * Math.exp(-rate * x);
		return density;
	}

}




class Gamma extends Distribution {

	initialise() {
		this.alpha = new Parameter(this.distNum, "alpha", "&alpha;", 1, 0, null, 0.2);
		this.beta = new Parameter(this.distNum, "beta", "&beta;", 1, 0, null, 0.2);
		this.parameters = [this.alpha, this.beta];
		this.norm = 1;
		this.desc = "Gamma distribution with shape alpha and scale beta. The mean is alpha/beta.";
	}

	getXRange() {
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		return [0, alpha/beta * 5];
	}

	prepareDensity() {
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		this.norm = beta ** alpha / math.gamma(alpha);
	}

	getDensity(x) {
		if (x <= 0) return 0;
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		var density = this.norm * Math.exp((alpha - 1) * Math.log(x) - (beta*x) );
		return density;
	}

}




class InverseGamma extends Distribution {

	initialise() {
		this.alpha = new Parameter(this.distNum, "alpha", "&alpha;", 1, 0, null, 0.2);
		this.beta = new Parameter(this.distNum, "beta", "&beta;", 1, 0, null, 0.2);
		this.parameters = [this.alpha, this.beta];
		this.norm = 1;
		this.desc = "Inverse-gamma distribution with shape alpha and scale beta.";
	}

	getXRange() {
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		return [0, alpha > 1 ? 5*beta/(alpha-1) : 15*beta/(alpha+1)];
	}

	prepareDensity() {
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		this.norm = beta ** alpha / math.gamma(alpha);
	}

	getDensity(x) {
		if (x <= 0) return 0;
		var alpha = this.alpha.get();
		var beta = this.beta.get();
		var density = this.norm * Math.exp((-alpha - 1) * Math.log(x) + (-beta/x) );
		return density;
	}

}




class Laplace extends Distribution {

	initialise() {
		this.rate = new Parameter(this.distNum, "rate", "&lambda;", 1, 0, null, 0.2);
		this.parameters = [this.rate];
		this.desc = "Laplace distribution with a rate of lambda, or a scale of 1/lambda.";
	}


	getXRange() {
		var rate = this.rate.get();
		return [-Math.log(2) / rate * 6, Math.log(2) / rate * 6];
	}

	getDensity(x) {
		var rate = this.rate.get();
		var density = rate/2 * Math.exp(-Math.abs(x) * rate);
		return density;
	}

}






