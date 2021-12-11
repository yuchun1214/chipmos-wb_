#include <include/route.h>
#include <cstdarg>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>

route_t::route_t()
{
    _wb_stations = {WB1, WB2, WB3, WB4};
    _da_stations = {DA1, DA2,  DA3,  DA4,  DA5,  DA6,  DA7, DA8,
                    DA9, DA10, DA11, DA12, DA13, DA14, DA15};
}

void route_t::setRoute(csv_t all_routes)
{
    // first, get the set of route names
    std::vector<std::string> route_names = all_routes.getColumn("route");
    std::set<std::string> route_list_set(route_names.begin(),
                                         route_names.end());
    route_names =
        std::vector<std::string>(route_list_set.begin(), route_list_set.end());

    // second, call route_t::setRoute for each route respectively
    csv_t df;
    iter(route_names, i)
    {
        df = all_routes.filter("route", route_names[i]);
        setRoute(route_names[i], df);
    }
}

void route_t::setRoute(std::string routename, csv_t dataframe)
{
    std::vector<std::vector<std::string> > data = dataframe.getData();
    std::vector<station_t> stations;
    std::map<std::string, std::string> elements;
    unsigned int size = dataframe.nrows();
    for (unsigned int i = 0; i < size; ++i) {
        elements = dataframe.getElements(i);
        stations.push_back(station_t{.route_name = elements["route"],
                                     .station_name = elements["desc"],
                                     .oper = stoi(elements["oper"]),
                                     .seq = stoi(elements["seq"])});
    }

    _routes[routename] = stations;

    setupBeforeStation(routename, true, 7, 4, WB1, WB2, WB3, WB4);  // WB - 7
}


void route_t::setQueueTime(csv_t queue_time_df)
{
    unsigned int nrows = queue_time_df.nrows();
    std::map<std::string, std::string> elements;
    std::map<int, int> queue_time;
    int station;
    for (unsigned int i = 0; i < nrows; ++i) {
        elements = queue_time_df.getElements(i);
        station = std::stoi(elements["station"]);
        elements.erase(elements.find("station"));
        for (std::map<std::string, std::string>::iterator it = elements.begin();
             it != elements.end(); it++) {
            queue_time[std::stoi(it->first)] = std::stoi(it->second) * 60;
        }
        _queue_time[station] = queue_time;
    }
}

std::vector<station_t> route_t::setupBeforeStation(std::string routename,
                                                   bool remove,
                                                   int nstations,
                                                   int nopts,
                                                   ...)
{
    std::set<int> opers;
    std::set<int> station_opers;
    std::set<unsigned int> indexes;

    std::vector<station_t> stations;

    va_list variables;
    va_start(variables, nopts);
    for (int i = 0; i < nopts; ++i) {
        opers.insert(va_arg(variables, int));
    }
    va_end(variables);


    int idx;
    for (unsigned int i = 0; i < _routes[routename].size(); ++i) {
        if (opers.count(_routes[routename][i].oper) !=
            0) {  // _routes[routename][i].oper is WB or DA
            idx = i - nstations;
            if (idx < 0) {
                idx = 0;
            }
            for (unsigned int j = idx; j <= i; ++j) {
                station_opers.insert(_routes[routename][j].oper);
                indexes.insert(j);
                // stations.push_back(_routes[routename][j]);
            }
        }
    }

    _beforeWB[routename] = station_opers;


    if (remove) {
        std::vector<int> vec(indexes.begin(), indexes.end());
        iter(vec, i) { stations.push_back(_routes[routename][vec[i]]); }
        _routes[routename] = stations;
    }

    return stations;
}

bool route_t::isLotInStations(lot_t lot)
{
    int idx, oper;
    if (_wb_stations.count(lot.oper()) &&
        lot.mvin()) {  // if lot is on WB station  and has moved in
        idx = findStationIdx(lot.route(),
                             lot.oper());  // locate the oper on the route
        if (idx > 0 && (unsigned int) (idx + 1) <
                           _routes[lot.route()]
                               .size()) {  // check if idx is route is resonable
            oper = _routes[lot.route()][++idx].oper;
            return _beforeWB[lot.route()].count(oper);
        } else
            return false;
    }
    return _beforeWB[lot.route()].count(lot.oper());
}

int route_t::findStationIdx(std::string routename, int oper)
{
    int idx = -1;
    iter(_routes[routename], i)
    {
        if (_routes[routename][i].oper == oper) {
            idx = i;
        }
    }

    return idx;
}


/**
 * In this function, the main issue is to determine the begining station of
 * traversal.
 *
 * if lot is in DA and lot isn't moved in, -> dispatch
 * if lot isn't in DA, ->
 *  1. traverse the route to DA and sum the queue time.
 *  2. dispatch
 * if lot is in DA and has been moved in ->
 *  1. begining station is next station.
 *
 */
int route_t::calculateQueueTime(lot_t &lot)
{
    if (lot.isTraversalFinish())
        return 0;

    std::string routename = lot.route();

    // locate the lot in the route
    int idx = findStationIdx(routename, lot.tmp_oper);

    // if can't locate the oper, idx will be less then 0
    if (idx < 0) {
        std::string error = "Can't locate lot.temp_oper(" +
                            std::to_string(lot.tmp_oper) + ") on the route";
        throw std::logic_error(error);
    }

    // check if lot is in WB
    // the conditions are as follow:
    //  1. if lot.oper == WB
    //  2. if lot.mvin == true
    //
    // if the conditions are true it means that lot is in WB and lot has moved
    // in ! in condition 2,  the lot is considered as lot is finished in WB. To
    // calculate the arrival time, lot is only need to plus the out plan time
    // and the queue time in next stations.
    //
    // in this condition, idx need to plus 1, the lot start from next station
    if (_wb_stations.count(lot.tmp_oper) == 1) {
        if (lot.tmp_mvin) {        // if lot is in WB an has moved in,
            idx += 1;              // move to next station
            lot.tmp_mvin = false;  // mvin = false
        } else {                   // lot is waiting at WB station
            lot.setTraverseFinished();
            return 0;
        }
    }

    // check if lot is in DA
    // if lot is in DA,
    if (_da_stations.count(lot.tmp_oper) == 1) {  // lot is in DA
        if (lot.tmp_mvin) {                       // lot is in DA and mvin
            idx += 1;                             // advance
            lot.tmp_mvin = false;
        } else {  // lot is in D/A and hasn't moved in
            lot.tmp_oper = _routes[routename][++idx].oper;  // advance
            return 2;  // advance and dispatch
        }
    }

    lot.tmp_mvin = false;


    // traverse the route from the begining station idx to D/A or W/B or to the
    // end of route if traverse to the end of route, it means that the lot will
    // not being scheduled to process in W/B station -> it is totally wrong!
    // please check route list or the route.h define macro. maybe _wb_stations
    // doesn't contain the extra W/B stations
    int oper;
    int prev = 0;
    int qt = 0;
    int times_of_passing_cure = 0;
    iter_range(_routes[routename], i, idx, _routes[routename].size())
    {
        oper = _routes[routename][i].oper;
        if (_queue_time.count(oper)) {  // check if oper is a big station?
            if (prev) {                 // prev is not null
                if (_queue_time[prev][oper] > 0) {
                    qt += _queue_time[prev][oper];
                    prev = oper;
                } else {
                    std::string error_text =
                        std::to_string(prev) + " -> " + std::to_string(oper) +
                        " is invalid queue time combination, please check "
                        "queue_time's input file.  ";

                    throw std::logic_error(error_text);
                }
            } else {
                prev = oper;
            }

            if (oper == CURE) {
                ++times_of_passing_cure;
            } else if (_da_stations.count(
                           oper)) {  // oper is a D/A station,  dispatch
                lot.tmp_mvin = false;
                lot.tmp_oper = _routes[routename][i + 1]  // ad
                                   .oper;  // lot traverse to DA station
                lot.addQueueTime(qt);
                return 1;
            } else if (_wb_stations.count(oper)) {  // traverse to W/B station
                lot.tmp_mvin = false;
                lot.tmp_oper = oper;
                lot.addQueueTime(qt);
                lot.setTraverseFinished();
                return 0;
            }
        }
    }

    return -1;
}
