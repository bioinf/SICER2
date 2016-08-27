
class Island_threshold :
	def __init__(self, total_tags, windowSize, gapSize, window_pvalue, genomeLength, bin_size, e_value_threshold):
		self.gap_size = gapSize/windowSize;
		self.genome_length = int(ceil(float(genomeLength)/windowSize));
		self.average = (total_tags * 1.0 / genomeLength) * windowSize;
		self.bin_size = bin_size;
		self.max_index = max (500, int(2*self.average));
		self.poisson_value = [];
		self.window_score = [];
		for index in xrange(self.max_index):
			prob = scipy.stats.poisson.pmf(index, self.average);
			self.poisson_value.append(prob);
			if ( index < self.average):
				self.window_score.append(0);
			else:
				self.window_score.append(-log(prob) if prob > 0 else 1000);
		self.max_index = len(self.poisson_value);
		self.min_tags_in_window = 0;
		sf = 1 ;
		while (sf > window_pvalue ):
			sf -= self.poisson_value[self.min_tags_in_window]
			self.min_tags_in_window += 1;
		self.gap_contribution = self.gap_factor();
		self.boundary_contribution = self.boundary();
		self.cumulative=[];
		prob = self.boundary_contribution * self.poisson_value[self.min_tags_in_window];
		score = -log(self.poisson_value[self.min_tags_in_window]);
		scaled_score = int(round(score/self.bin_size));
		self.island_expectation =[0] * (scaled_score+1);
		self.island_expectation[scaled_score] = prob*self.genome_length;
		self.island_expectation[0] = self.boundary_contribution*self.genome_length/self.gap_contribution;
		self.threshold = self.find_island_threshold(e_value_threshold)

	def single_gap_factor(self):
		my_gap_factor=0;
		for i in xrange(self.min_tags_in_window):
			my_gap_factor +=self.poisson_value[i];
		return my_gap_factor;

	def gap_factor(self):
		if self.gap_size == 0 : return 1
		i = 1;
		gap_contribution = 1;
		my_gap_factor = self.single_gap_factor();
		for i in range(1, self.gap_size+1): gap_contribution += pow(my_gap_factor, i);
		return gap_contribution;

	def boundary(self):
		temp = self.single_gap_factor();
		temp = pow(temp, self.gap_size+1); 
		return temp*temp; # start & end 

	def background_island_expectation (self, scaled_score):
		current_max_scaled_score = len(self.island_expectation)-1;
		if scaled_score > current_max_scaled_score:
			#index is the scaled_score
			for index in range(current_max_scaled_score + 1, scaled_score+1):
				temp=0.0;
				#i is the number of tags in the added window
				i = self.min_tags_in_window;
				while ( int(round(index - self.window_score[i]/self.bin_size))>=0):
				#while ( (index - self.window_scaled_score[i])>=0):
					temp += self.poisson_value[i]* self.island_expectation[int(round(index - self.window_score[i]/self.bin_size))];
					#temp += self.poisson_value[i]* self.island_expectation[index - self.window_scaled_score[i]];
					i += 1;
				temp *= self.gap_contribution;
				self.island_expectation.append(temp);
				#print index, temp, self.island_expectation[index];
		return self.island_expectation[scaled_score];

	def find_island_threshold(self, e_value_threshold):
		threshold = .0000001*e_value_threshold;
		current_scaled_score = len (self.island_expectation) - 1;
		current_expectation = self.island_expectation[-1];
		assert (current_expectation == self.island_expectation[current_scaled_score]);
		interval = int(1/self.bin_size);
		if len(self.island_expectation) > interval:
			partial_cumu = sum(self.island_expectation[-interval: -1])
		else:
			partial_cumu = sum(self.island_expectation)
		while ( partial_cumu > threshold or  partial_cumu <1e-100):
			current_scaled_score += interval;
			current_expectation=self.background_island_expectation(current_scaled_score);
			if len(self.island_expectation) > interval:
				partial_cumu = sum(self.island_expectation[-interval: -1])
			else:
				partial_cumu = sum(self.island_expectation)

		# generate_cumulative_dist
		self.cumulative=[0]*len(self.island_expectation);
		partial_sum = 0.0
		for index in range(1, len(self.island_expectation)+1):
			complimentary = len(self.island_expectation) - index;
			partial_sum += self.island_expectation[complimentary]; # The end is outside of the index
			self.cumulative[complimentary]=partial_sum;

		for index in xrange(len(self.cumulative)):
			if self.cumulative[index]<=e_value_threshold:
				score_threshold = index*self.bin_size;
				break;
		return score_threshold;
