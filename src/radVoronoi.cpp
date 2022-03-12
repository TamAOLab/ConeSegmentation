#include "radVoronoi.h"

#include <float.h>
#define JC_VORONOI_IMPLEMENTATION
#define JCV_REAL_TYPE double
#define JCV_ATAN2 atan2
#define JCV_SQRT sqrt
#define JCV_FLT_MAX DBL_MAX
#define JCV_PI 3.141592653589793115997963468544185161590576171875
//define JCV_EDGE_INTERSECT_THRESHOLD 1.0e-10F
#include "jc_voronoi.h"

static int opcode(const jcv_rect& bb, const jcv_point& pt)
{
	int op = 0;
	if (pt.x <= bb.min.x) op |= 1;
	if (pt.x >= bb.max.x) op |= 2;
	if (pt.y <= bb.min.y) op |= 4;
	if (pt.y >= bb.max.y) op |= 8;
	return op;
}

std::vector<DoublePointArray> VoronoiEdges(DoublePointArray& markers, int width, int height)
{
	std::vector<DoublePointArray> res;

	if (markers.size() < 3) return res;

	double xmax = double(width) - 1.001;
	double ymax = double(height) - 1.001;
	jcv_rect bounding_box = { { 0.001, 0.001 }, { xmax, ymax } };

	jcv_point* points = new jcv_point[markers.size()];
	for (int i = 0; size_t(i) < markers.size(); i++) {
		DoublePointType& pt = markers[i];
		points[i].x = pt[0];
		points[i].y = pt[1];
	}

	jcv_diagram diagram;
	memset(&diagram, 0, sizeof(jcv_diagram));
	jcv_diagram_generate(int(markers.size()), (const jcv_point*)points, &bounding_box, NULL, &diagram);

	const jcv_edge* edge = jcv_diagram_get_edges(&diagram);
	while (edge)
	{
		if ((opcode(bounding_box, edge->pos[0]) & opcode(bounding_box, edge->pos[1])) == 0)
		{
			res.resize(res.size() + 1);
			DoublePointArray& pts = res[res.size() - 1];
			pts.resize(2);
			pts[0][0] = edge->pos[0].x;
			pts[0][1] = edge->pos[0].y;
			pts[1][0] = edge->pos[1].x;
			pts[1][1] = edge->pos[1].y;
		}
		edge = jcv_diagram_get_next_edge(edge);
	}

	delete[] points;

	return res;
}

