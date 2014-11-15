#ifndef __AVY__DEBUG__HPP_
#define __AVY__DEBUG__HPP_
#include <string>
#include <set>

namespace avy
{

#ifndef NAVYLOG
#define LOG(TAG,CODE) do { if (::avy::AvyLogFlag && ::avy::AvyLog.count (TAG) > 0) { CODE; } } while (0)
extern bool AvyLogFlag;
extern std::set<std::string> AvyLog;
void AvyEnableLog (std::string x);
#else
void AvyEnableLog (std::string x);
#define LOG(TAG,CODE) do { } while (0)
#endif

}


#endif
