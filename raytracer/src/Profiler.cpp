#include "Profiler.h"

map<string, Timer> Profiler::idToTimer = map<string, Timer>();
map<string, long> Profiler::idToCounter = map<string, long>();

ostream& operator<<(ostream& os, const Timer t)
{
	os << "Time: " << t.time <<"s, count: " << t.count << ", average time: " << (t.time / t.count);
	return os;
};

void Profiler::printReport()
{
	cout << "Profiling report:" << endl;
	for(map<string, Timer>::iterator iter = idToTimer.begin(); iter != idToTimer.end(); ++iter)
	{
		cout << "---" << iter->first << "::::" << iter->second << endl;
	}
	cout << "------------------------------------------" << endl;
	for(map<string, long>::iterator iter = idToCounter.begin(); iter != idToCounter.end(); ++iter)
	{
		cout << "---" << iter->first << "::::" << iter->second << endl;
	}
	cout << "End of profiling report!!!" << endl;
}

void Profiler::finishTimer( string id )
{
	idToTimer[id].timer.endTimer();
	idToTimer[id].time += idToTimer[id].timer.elapsedTime();
	idToTimer[id].count ++;
}

void Profiler::startTimer( string id )
{
	//iii ++;
	idToTimer[id].timer.beginTimer();
}


void Profiler::increaseCounter(string id, long i)
{
	idToCounter[id] += i;
}


void Profiler::clear()
{
	idToTimer.clear();
	idToCounter.clear();
}

#ifdef _DEBUG

#else
	// make sure that all profilers deleted from other places	
	;lkajsdf;j;
#endif
