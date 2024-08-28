var radius = function( d) { 

	if( d.radius )
		return d.radius;

	var radius;
	switch(d.group) {
    case 1:
        radius = 100;
        break;
    case 2:
        radius = 40;
        break;
    case 3:
        radius = 30;
        break;
    default:
        radius = 30;
				break;
	};

	return radius;
};

function replaceContentInContainer(target, source) {
	document.getElementById(target).innerHTML =  document.getElementById(source).innerHTML;
}


var menu_tree = function(selector, w, h) {

	// Clear all svg
	d3.select(selector).selectAll('svg').remove();

	// Create the svg 
  this.svg = d3.select(selector).append('svg:svg')
					.attr('width', w )
					.attr('height', h );

	// Get bounding box measures for the svg
	var outer_bbox = this.svg.node().getBoundingClientRect();

	this.outer_w = outer_bbox.width;
	this.outer_h = outer_bbox.height;

	// Create the rectangle canvas
	this.svg.append('svg:rect')
		.style('stroke', '#fff')
		.style('fill', '#fff')
		.attr('opacity',0.2)
	  .attr('x', this.outer_w*0.025)
	  .attr('y', this.outer_h*0.025)
	  .attr('width', this.outer_w*0.95)
	  .attr('height', this.outer_h*0.95)
		.attr('id','canvas');

	// Get the bounding box for the inner rectangle object
	var box = d3.select('#canvas').node().getBoundingClientRect();
		
	this.w = outer_bbox.width;
	this.h = outer_bbox.height;
	
	// Create the forced tree layout
  this.force = d3.layout.force()
    .on('tick', this.tick.bind(this))
    .charge(function(d) { return -10*radius(d); })
    .linkDistance(function(d) { return d.source.distance ?  150 : ( d.target.group == 3 ? 200 : 200) ; })
    .size([ this.w, this.h]);

};

menu_tree.prototype.update = function(json) {

  if (json) this.json = json;

  this.json.fixed = true;
  this.json.x = this.w/2;
  this.json.y = this.h/2; 

	var nodes;
	var links;

	if(0){
	  nodes = this.json.nodes;
		links = this.json.links;
	} else {
		nodes = this.flatten(this.json);
	  links = d3.layout.tree().links(nodes);
	}
  var total = nodes.length || 1;
	this.nodes = nodes;

  // Restart the force layout
  this.force
    .nodes(nodes)
    .links(links)
    .start();

  // Update the links
  this.link = this.svg.selectAll('line.link')
    .data(links, function(d) { return d.target.name; });

  // Enter any new links
  this.link.enter().insert('svg:line', '.node')
    .attr('class', 'link')
    .attr('x1', function(d) { return d.source.x; })
    .attr('y1', function(d) { return d.source.y; })
    .attr('x2', function(d) { return d.target.x; })
    .attr('y2', function(d) { return d.target.y; })
		.attr('color','rgba(255,255,255,0.5)');

  // Exit any old links.
  this.link.exit().remove();

	// Update the nodes
  this.node = this.svg.selectAll('.node')
    .data(nodes, function(d) { return d.name; });
	
	// Enter any new nodes
  this.eltEnter = this.node.enter();
 
  this.eltEnter.append('svg:circle')
			.attr('class', 'node')
	    .attr('r', function(d) { return radius(d); })
		  .style('fill',  'rgba(245,128,37,1)' )
			.call(this.force.drag)
			.on('click', this.click.bind(this));
 
	this.eltEnter.append('foreignObject')
		.attr('class', 'node')
		.attr('width', function(d) {return 1.75*radius(d) ; })
    .attr('height', function(d) {return radius(d) ; })
		.attr('x',function(d) {return -radius(d)*1.75/2 ; })
		.attr('y',function(d) {return -radius(d)/4 ; })
		.append('xhtml:div')
		.attr('class','node')
		.attr('display', 'table')
		.attr('width', function(d) {return 1.75*radius(d) ; })
    .attr('height',function(d) {return 2*radius(d) ; })
    .append('xhtml:p')
		.attr('class', 'node')
		.attr('display','table-cell')
		.attr('text-align','center')
		.attr('vertical-align','middle')
		.attr('align','center')
		.style('text-anchor','middle')
		.style('color','white')
    .style('font-family','Source Sans Pro')
		.style('font-size', function(d) { return d.fontsize ? d.fontsize : Math.min(0.5* radius(d)*0.925,100)*2/3 + 'px'; })
		.style('line-height',function(d) { return Math.min(0.5* radius(d),100)/2 + 'px'; })
    .html(function(d){ return d.name ;} )
		.call(this.force.drag)	;
 
  this.node.attr('fixed',true);
     
  return this;
};

menu_tree.prototype.flatten = function(root) {
  var nodes = [], i = 0;

  function recurse(node) {
    if (node.children) {
      node.size = node.children.reduce(function(p, v) {
        return p + recurse(v);
      }, 0);
    }
    if (!node.id) node.id = ++i;
    nodes.push(node);
    return node.size;
  }

  root.size = recurse(root);
  return nodes;
};

menu_tree.prototype.click = function(d) {
	
	//replaceContentInContainer('profile','empty_profile');
	if( d.divid == "home"  ) {

		d3.select('#welcome').style('display','block');

		d3.select('#menu_tree_main').style('display','block');

		d3.select('#menu_tree_side').style('display','none');

		replaceContentInContainer('main','empty_profile');

		d3.select('#profile').style("display", "block");

		
	} else {
		
		if (d.divid) {
		replaceContentInContainer('main',d.divid);
		} else {
		replaceContentInContainer('main','publications');
		};

		d3.select('#menu_tree_main').style('display','none');
		d3.select('#menu_tree_side').style('display','block');
		d3.select('#welcome').style('display', 'none');
		d3.select('#profile').style('display', 'none');

		//d3.select('#menu_bar').transition()
		//  .style("height", "auto");
	};

};

menu_tree.prototype.mouseover = function(d) {

};

menu_tree.prototype.mouseout = function(d) {
 
};

menu_tree.prototype.tick = function() {
  var h = this.h;
  var w = this.w;

	var q = d3.geom.quadtree(this.nodes),
      i = 0,
      n = this.nodes.length;

  while (++i < n) q.visit(collide(this.nodes[i]));

	this.node.attr('cx', function(d) { return d.x =  Math.max(radius(d) , Math.min(w - radius(d) , d.x)); })
        .attr('cy', function(d) { return d.y =  Math.max(radius(d), Math.min(h - radius(d), d.y)); });
	
  this.svg.selectAll('.node')
		.attr('transform', function(d) {
    return 'translate(' + Math.max(radius(d), Math.min(w - radius(d), d.x)) + ',' +
		 	Math.max(radius(d), Math.min(h - radius(d), d.y)) + ')'; });

  this.link.attr('x1', function(d) { return d.source.x; })
    .attr('y1', function(d) { return d.source.y; })
    .attr('x2', function(d) { return d.target.x; })
    .attr('y2', function(d) { return d.target.y; });
	return this;
};

menu_tree.prototype.cleanup = function() {
  this.update([]);
  this.force.stop();
};


function collide(node) {
var r = node.radius + 16,
    nx1 = node.x - r,
    nx2 = node.x + r,
    ny1 = node.y - r,
    ny2 = node.y + r;
	return function(quad, x1, y1, x2, y2) {
    if (quad.point && (quad.point !== node)) {
        var x = node.x - quad.point.x,
            y = node.y - quad.point.y,
            l = Math.sqrt(x * x + y * y),
            r = node.radius + quad.point.radius;
        if (l < r) {
            l = (l - r) / l * .5;
            node.x -= x *= l;
            node.y -= y *= l;
            quad.point.x += x;
            quad.point.y += y;
        }
    }
    return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
	};
}




