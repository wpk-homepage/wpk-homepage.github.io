'use strict';

function plotdata(data) {
	nv.addGraph(function() {
		var chart = nv.models.scatterChart()
			//.showDistX(true)
			//.showDistY(true)
			.color(d3.scale.category10().range());

		chart.xAxis.tickFormat(d3.format('2d')).axisLabel('Helling (%)');
		chart.yAxis.tickFormat(d3.format('2d')).axisLabel('Snelheid (km/u)');
		
		chart.xDomain([-25, 25]);
		chart.yDomain([0, 40]);

		d3.select('#chart svg')
			.datum(data)
			.transition().duration(500)
			.call(chart);

		nv.utils.windowResize(chart.update);

		return chart;
	});
}

