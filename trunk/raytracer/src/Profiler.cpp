#include "Profiler.h"

map<string, Timer> Profiler::idToTimer = map<string, Timer>();

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

void Profiler::clear()
{
	idToTimer.clear();
}

#ifdef _DEBUG

#else
	// make sure that all profilers deleted from other places	
	;lkajsdf;j;
#endif
