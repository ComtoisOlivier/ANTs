/*
 * ANTSUseLandmarkImagesWithTextFileToGetThinPlateDisplacementField.h
 *
 *  Created on: Jul 20, 2015
 *      Author: olcoma
 */

#ifndef EXAMPLES_INCLUDE_ANTSUSELANDMARKIMAGESWITHTEXTFILETOGETTHINPLATEDISPLACEMENTFIELD_H_
#define EXAMPLES_INCLUDE_ANTSUSELANDMARKIMAGESWITHTEXTFILETOGETTHINPLATEDISPLACEMENTFIELD_H_

namespace ants
{
extern int ANTSUseLandmarkImagesWithTextFileToGetThinPlateDisplacementField( std::vector<std::string>, // equivalent to argv of command line
                                                                                  // parameters to main()
                                                        std::ostream* out_stream  // [optional] output stream to write
                                                        );
} // namespace ants

#endif /* EXAMPLES_INCLUDE_ANTSUSELANDMARKIMAGESWITHTEXTFILETOGETTHINPLATEDISPLACEMENTFIELD_H_ */
