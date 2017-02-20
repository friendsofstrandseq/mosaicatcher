#include <iostream>
#include <vector>
#include "hmm.h"
#include "seqTracks.h"
#include "trellis.h"

/* Notes 1:
	If I decide to use GHMM, check out this post:
		https://sourceforge.net/p/ghmm/mailman/message/822875/
*/

int main()
{
	StochHMM::model hmm;
	StochHMM::StateFuncs default_functions;
	std::string s = "src/strseq.gaussian.hmm";
	hmm.import(s,&default_functions);
	std::cout << "successfully loaded model " << s << std::endl;
	hmm.print();

	// Input (one strand only for now)
	std::vector<double> plus = {0.1225978, 0.5659937, -0.1685184, -0.1748695, -0.04298955, 0.1899047, 0.01009715, -0.02918318, 0.1902079, -0.06849704, 0.15036, -0.04167554, 0.003029844, -0.002109414, 0.1842016, -0.2634236, -0.3295542, 0.305007, 0.04624348, 0.03049509, 0.9272208, 1.277645, 0.6107512, 1.171286, 1.174998, 0.9974216, 0.8529463, 1.191086, 0.9969731, 0.8537285, 0.9963578, 1.269444, 1.073113, 1.259121, 0.8568738, 0.9288473, 1.046971, 1.273443, 1.015991, 1.138565, 0.8396678, 0.8172367, 0.7810229, 1.146505, 0.8251086, 0.9141333, 0.9599123, 0.6829162, 1.256244, 1.097067, 0.3756536, 0.6862918, 0.278008, 0.309078, 0.579166, 0.5184487, 0.4589595, 0.2589538, 0.3186686, 0.2198736, 0.4460325, 0.5318255, 0.4446543, 0.4587801, 0.5757608, 0.5691747, 0.5253653, 0.5679433, 0.6843114, 0.3285339, 0.3296672, 0.29591, 0.4074933, 0.3650981, 0.7671388, 0.3877953, 0.4933179, 0.6765814, 0.555503, 0.6661862, 0.5845685, 0.4026926, 0.3292454, 0.61933, 0.1761858, 0.3795035, 0.6490185, 0.4803516, 0.06274393, 0.3514703, 0.6008354, 0.4855625, 0.4349988, 0.6140038, 0.5600466, 0.3775088, 0.3103108, 0.45385, 0.4301768, 0.3684095};

	// Input sequence into HMM
	StochHMM::tracks* trks = hmm.getTracks();
	StochHMM::track*  trk  = trks->getTrack("PLUS");
	StochHMM::sequences sqs(trks->size());
	StochHMM::sequence* sq = new(std::nothrow) StochHMM::sequence(&plus,trk);
	sq->print();

	// Set up a "trellis", whatever that means
	StochHMM::trellis trell(&hmm, &sqs);
	trell.viterbi();
	// >< >< crash >< ><

}
