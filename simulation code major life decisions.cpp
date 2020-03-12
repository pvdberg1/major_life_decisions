#include <iostream>
#include <random>
#include <ctime>
#include <fstream>
#include <vector>

using namespace std;

/////////////
// GLOBALS //
/////////////

struct option {
double obj;
double subj;
int eval;
};

vector<option> consider;				// the options currently under consideration
vector<int> evalvec;					// this vector stores, per parameter combination, for how many timesteps each option is exploited

// parameter combinations for treshold 1 & 2, the accuracy, and the sample probability:
const double minthresh1 = 0.0;
const double stepthresh1 = 1.0;
const int numstepthresh1 = 2;

const double minthresh2 = 1.0;
const double stepthresh2 = 0.05;
const int numstepthresh2 = 31;

const double minnoise = 0.5;
const double stepnoise = 1.0;
const int numstepnoise = 2;

const double minprobsample = 0.0;
const double stepprobsample = 0.01;
const double numstepprobsample = 21;

const int totalcombs = 2604;			//total number of parameter combinations (compute this yourself)

//fixed parameters
double chancetocheck = 0.5;             //chance to check option when about to make choice
const double recency = 0.5;				//the relative weight of the past versus the present
const int timelimit = 100;				//total number of time steps available
int const maxsize = 3;					//maximum number of options to consider at one time

const int reps = 1000;					//number of replicates

double data1[12][totalcombs];			//datavector that is eventually outputted
//vector<vector<double>> bigdata;		//to output the evaluation time for each run

//random stuff
mt19937 mt(time(0));
uniform_int_distribution<int> Random(0, (maxsize-1));
uniform_real_distribution<double> Uniform(0, 1);
normal_distribution<double> Normal(0, 1);

ofstream output;						//output over all parameter combinations
ofstream output2;						//output of single run
ofstream output_eval;					//output of all options and the number of times they have been evaluated

int cnt_reps = 0;

//////////////
// FUNTIONS //
//////////////

int explore(double noise, double thresh1, double thresh2, double probsample)
{
	
	option thisopt;
    //the objective value of a new option is a random number from a normal distribution
    thisopt.obj = Normal(mt);

	//the subjective value of a new option is equal to the objective value and some noise
	thisopt.subj = thisopt.obj + Normal(mt)*noise;
    thisopt.eval = 1;

    
	//the objective value is only added to S if it exceeds thresh1
	if(thisopt.subj > thresh1)
	{
		int minsubj = 0;
			
		//it is added to S if S does not yet have maximum size...
		if(consider.size() < maxsize)
		{
			consider.push_back(thisopt);
		}
		//or if the option in S with the smallest subjective value has smaller subjective value than the new option

		
		else
		{
			for(int i = 1; i < consider.size(); ++i)
			{
				if(consider[i].subj <= consider[minsubj].subj) {minsubj = i;}
			}

			if(consider[minsubj].subj < thisopt.subj)
			{
				//vector<double> thisvec;

				//thisvec.push_back(cnt_reps);
				//thisvec.push_back(consider[minsubj].eval);
				//thisvec.push_back(thresh1);
				//thisvec.push_back(thresh2);
				//thisvec.push_back(noise);
				//thisvec.push_back(probsample);

				//bigdata.push_back(thisvec);

				evalvec.push_back(consider[minsubj].eval);

				//thisvec.clear();

				consider[minsubj] = thisopt;
			}
		}
		//if the new option has a higher subjective value than tresh 2, this function returns its index as the chosen option
		if(thisopt.subj > thresh2) {return minsubj;}
		//if this is not the case, this function returns -1, which means 'no choice made'
		else {return -1;}
	}
	
	return -1;
}

//int exploit(double thresh1, double noise, int whathappened) <- only used for sims to check dynamics of single runs
int exploit(int choicemade, double noise, double thresh1, double thresh2)
{
    int pick;

    //if no choice is made, pick a random option in S for exploitation
    if (choicemade == -1)
    {
        //draw a random integer between 0 and the size of S, to draw a random option
        uniform_int_distribution<int> RandomOpt(0, consider.size() - 1);
        pick = RandomOpt(mt);
    }

    //otherwise exploit the option that is about to be chosen
    else
    {
        pick = choicemade;
    }

	//update the option according to exploitation process
    consider[pick].subj = (consider[pick].subj * consider[pick].eval * recency + consider[pick].obj + Normal(mt) * noise)/(consider[pick].eval * recency + 1);
    consider[pick].eval += 1;

	//if the subjective value of the newly exploited option came below tresh 1, remove it from S
	if(consider[pick].subj < thresh1)
    {
        consider.erase(consider.begin()+pick);
    }

	//if the exploited option has a higher subjective value than tresh 2, this option's index becomes the choicemade
	if(consider[pick].subj > thresh2)
    {
        choicemade = pick;
    }
	//otherwise, choicemade becomes -1 (= no choice made)
	else
    {
        choicemade = -1;
    }

	return(choicemade);

	//this is only used for sims to check dynamics of single runs:
	//whathappened = pick + 1;
	//return whathappened;
}

//this calculates averages and variances over all replicates of a single parameter combination
void stats(double data2[3][reps], int cnt)
{
    for(int i = 0; i<3; ++i)
    {
        double sum = 0;
        double avg = 0;
        double ss = 0;
        double var = 0;

        for(int j = 0; j < reps; ++j)
        {
            sum += data2[i][j];
        }

        avg = sum/reps;

        for(int j = 0; j < reps; ++j)
        {
            ss += (avg - data2[i][j]) * (avg - data2[i][j]);
        }

        var = ss/reps;

        double sdev = sqrt(var);

        data1[i*2][cnt] = avg;
        data1[i*2+1][cnt] = sdev;
    }

	  //gem en sd van eval vector
	  double sum = 0;
      double avg = 0;
      double ss = 0;
      double var = 0;

      for(int j = 0; j < evalvec.size(); ++j)
      {
          sum += evalvec[j];
      }

      avg = sum/evalvec.size();

      for(int j = 0; j < evalvec.size(); ++j)
      {
          ss += (avg - evalvec[j]) * (avg - evalvec[j]);
      }

      var = ss/evalvec.size();

      double sdev = sqrt(var);

      data1[6][cnt] = avg;
      data1[7][cnt] = sdev;
}

//this is only used to output dynamics of single runs
void writeheaders2()
{
    output2 << "opt1_subj\t" << "opt1_obj\t" << "opt1_eval\t"
     << "opt2_subj\t" << "opt2_obj\t" << "opt2_eval\t"
     << "opt3_subj\t" << "opt3_obj\t" << "opt3_eval\t" << "treshold1\t" << "treshold2\t"<< "choicemade\t" << "whathappened\n";
}

// the data are written here
void writedata(int numcombs)
{
    output << "timestep_avg\t" << "timestep_sd\t" << "obj_avg\t" << "obj_sd\t" << "subj-obj_avg\t" << "subj-obj_sd\t"
     << "eval_avg\t"<< "eval_sd\t"<< "thresh1\t" << "thresh2\t" << "noise\t" << "probsample\n";

    for(int i = 0; i < numcombs; ++i)
    {
        for(int j = 0; j < 11; ++j)
        {
             output << data1[j][i] << "\t";
        }

        output << data1[11][i] << "\n";
    }
}

/*
void writedata_eval()
{
	 output_eval << "run\t" << "eval\t" << "thresh1\t" << "thresh2\t" << "noise\t" << "probsample\n";

	 for(int i = 0; i < bigdata.size(); ++i)
	 {
		 for(int j = 0; j < 5; ++j)
		 {
			output_eval << bigdata[i][j] << "\t";
		 }
		 output_eval << bigdata[i][5] << "\n";
	 }
}
*/

int main()
{
    int cnt = 0;
	
    for(int m = 0; m < numstepthresh1; ++m)
    {
        for(int n = 0; n < numstepthresh2; ++n)
        {
            for(int o = 0; o < numstepnoise; ++o)
            {
                for(int p = 0; p < numstepprobsample; ++p)
                {
					double thresh1 = minthresh1 + m * stepthresh1;
					double thresh2 = minthresh2 + n * stepthresh2;
					double noise = minnoise + o * stepnoise;
					double probsample = minprobsample + p * stepprobsample;
					
					//we only run this parameter combination if treshold 2 is larger than treshold 1
					if(thresh1 <= thresh2)
                    {
			            //this is a matrix to store the main outputs of all the replicates within this parameter combination
						double data2[3][reps];

						if(cnt % 100 == 0) {cout << cnt << "\n";}

                        for(int k = 0; k < reps; ++k)
                        {
							cnt_reps += 1;

							//this is only for outputting the dynamics of a single run:
							//output2.open("output_singlerun.txt");
							//writeheaders2();

							//follows the timestep we're in
							int step = 0;
							//checks if we have made the choice; -1 means 'no choice made'
							int choicemade = -1;
							//clear the vector of options
							consider.clear();

							// this runs for as long as a choice has not yet been made (later in this loop, a choice is forced)
							// if 'step' is equal to the maximum number of allowed timesteps.
							while(choicemade < 0)
                            {
                                ++step;

								//to follow what happens in each timestep of a single run. Either becomes equal to an index
								//of one of the options (meaning that that option has been exploited), or becomes equal to the
								//maximum size of S, meaning that exploration happened.
								//int whathappened = -1;

								//must explore if there are no options worth considering in S
								if(consider.size() == 0)
								{
									choicemade = explore(noise, thresh1, thresh2, probsample);
									//whathappened = maxsize;
								}
								else
								{
									//explore if a randomly drawn number is lower than parameter 'probsample'
									if(Uniform(mt) < probsample)
									{
										choicemade = explore(noise, thresh1, thresh2, probsample);
										//whathappened = maxsize;
									}
									//otherwise exploit
									else
									{
										choicemade = exploit(choicemade, noise, thresh1, thresh2);
										//whathappened = exploit(thresh1, noise, whathappened);
									}
								}

								//if a choice has been made, there is a chance to exploit chosen option depending on parameter 'chancetocheck'
								//in this case, 'choicemade' contains the index of the chosen option in S
								while(choicemade > -1)
								{
                                    if (Uniform(mt) < chancetocheck)
                                    {
                                        choicemade = exploit(choicemade, noise, thresh1, thresh2);
                                    }
                                    //if the chosen option is not further exploited, document in which timestep the choice was made, the subjective and objective value of chosen option, and terminate while loop
                                    else
                                    {
										data2[0][k] = step;
										data2[1][k] = consider[choicemade].obj;
										data2[2][k] = consider[choicemade].subj - consider[choicemade].obj;
										break;
                                    }
								}

								//if we got to the last timestep without making a choice, a choice is forced
								if(step == timelimit)
								{
									data2[0][k] = step;

									// in the version below, the option with the best perceived fit gets chosen if a choice is not made at the end.
									//if there are no options in S, a completely random choice is made.
									double maxsubj = -1000;
									int bestopt = -1;

									if(consider.size() > 0)
									{

                                        for(int j = 0; j < consider.size(); ++j)
                                        {
                                            if(consider[j].subj > maxsubj)
                                            {
                                                maxsubj = consider[j].subj;
                                                bestopt = j;
                                            }
                                        }

									data2[1][k] = consider[bestopt].obj;
									data2[2][k] = consider[bestopt].subj - consider[bestopt].obj;
									}

									else
									{
										data2[1][k] = Normal(mt);
										data2[2][k] = data2[1][k] + Normal(mt)*noise;
									}

									choicemade = 0;


									//in this version, a random choice between the options in S is made if no choice has been made at the end.
									//because there is no logical order to the choices, a random choice is just implemented as the first option in S.
									//if there are no options in S, a completely random choice is made.
									/*if(consider.size() > 0)
									{
										data2[1][k] = consider[0].obj;
										data2[2][k] = consider[0].subj - consider[0].obj;
									}
									else
									{
										data2[1][k] = Normal(mt);
										data2[2][k] = data2[1][k] + Normal(mt)*noise;
									}

									choicemade = 0;*/
								}

								//write time series data for a the very first run only
								//if(k == 0 && cnt == 0)
								//{
									//writedata2(thresh1, thresh2, whathappened, choicemade);
								//}

                            }

                            //output2.close();
							
							for (int i = 0; i < consider.size(); ++i)
							{
								//vector<double> thisvec;

								//thisvec.push_back(cnt_reps);
								//thisvec.push_back(consider[i].eval);
								//thisvec.push_back(thresh1);
								//thisvec.push_back(thresh2);
								//thisvec.push_back(noise);
								//thisvec.push_back(probsample);

								//bigdata.push_back(thisvec);

								evalvec.push_back(consider[i].eval);

								//thisvec.clear();
							}
	                    }

                        data1[8][cnt] = thresh1;
                        data1[9][cnt] = thresh2;
                        data1[10][cnt] = noise;
                        data1[11][cnt] = probsample;


                        // calculate averages and standard deviations of response variables over all replicates of this parameter combo
						stats(data2, cnt);
						cnt += 1;
                    }
					evalvec.clear();
				}
            }
        }
    }
	//write output to file
	output.open("output.txt");
	//output_eval.open("output_eval.txt");
    writedata(cnt);
	//writedata_eval();

    return 0;
}
