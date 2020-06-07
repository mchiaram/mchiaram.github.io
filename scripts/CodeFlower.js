var radius = function( d) { 
	var radius;

	switch(d.group) {
    case 1:
        radius = 150;
        break;
    case 2:
        radius = 50;
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


var CodeFlower = function(selector, w, h) {

  d3.select(selector).selectAll('svg').remove();

  this.svg = d3.select(selector).append('svg:svg')
				.attr('width', '100vw' )
    .attr('height', '100vh' )
		.style('fill', '#000');


  this.canvas = this.svg.append('svg:rect')
    .style('stroke', '#000')
    .style('fill', '#fff')
    .attr('width', '100vw' )
    .attr('height', '100vh' )
		.attr('id','canvas');

	this.canvas.append('text').attr('x',0).attr('y',0).text('Work in progress');

	var box = d3.select(selector).node().getBoundingClientRect();
	this.w = box.width;
	this.h = box.height;
	
  this.force = d3.layout.force()
    .on('tick', this.tick.bind(this))
    .charge(function(d) { return -50*radius(d); })
    .linkDistance(function(d) { return d.source.distance ?  d.source.distance : ( d.target.children ? 150 : 50) ; })
    .size([ this.w, this.h]);

};

CodeFlower.prototype.update = function(json) {

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

  // remove existing text (will readd it afterwards to be sure it's on top)
  this.svg.selectAll('text').remove();

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
  this.node = this.svg.selectAll('circle.node')
    .data(nodes, function(d) { return d.name; });
	
	// Enter any new nodes
  var eltEnter =  this.node.enter();
 
  var cicles = eltEnter.append('svg:circle')
			.attr('class', 'node')
	    .attr('r', function(d) { return radius(d); })
		  .style('fill', 'rgba(245,128,37,1)')
			.call(this.force.drag)
	    .on('click', this.click.bind(this))
		  .on('mouseover', this.mouseover.bind(this))
			.on('mouseout', this.mouseout.bind(this));
 
    eltEnter.append('svg:text')
      .attr('class', 'node')
			.text(function (d) { return d.name;})
			.attr('text-anchor', 'middle')
			.attr('font-family','Montserrat')
			.attr('fill','#fff')
			.attr('font-size',function (d) { return d.children ? '16' : '10'; } )
			.each(getSize)
			.call(this.force.drag)
			.on('click', this.click.bind(this))
			.on('mouseover', this.mouseover.bind(this))
			.on('mouseout', this.mouseout.bind(this));
 
    this.node.attr('fixed',true);
     
  // Exit any old nodes
  this.node.exit().remove();

  return this;
};

CodeFlower.prototype.flatten = function(root) {
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

CodeFlower.prototype.click = function(d) {
	window.location = d.url;
};

CodeFlower.prototype.mouseover = function(d) {


};

CodeFlower.prototype.mouseout = function(d) {
 // this.text.style('display', null);
};

CodeFlower.prototype.tick = function() {
  var h = this.h;
  var w = this.w;
	

	var q = d3.geom.quadtree(this.nodes),
      i = 0,
      n = this.nodes.length;

  while (++i < n) q.visit(collide(this.nodes[i]));


	this.node.attr('cx', function(d) { return d.x = Math.max(radius(d) , Math.min(w - radius(d) , d.x)); })
        .attr('cy', function(d) { return d.y = Math.max(radius(d), Math.min(h - radius(d), d.y)); });
	
  this.svg.selectAll('circle.node')
		.attr('transform', function(d) {
    return 'translate(' + Math.max(radius(d), Math.min(w - radius(d), d.x)) + ',' +
		 	Math.max(radius(d), Math.min(h - radius(d), d.y)) + ')';
		});

  this.svg.selectAll('text.node')
		.attr('transform', function(d) {
    return 'translate(' + Math.max(radius(d), Math.min(w - radius(d), d.x)) + ',' + Math.max(radius(d), Math.min(h - radius(d), d.y)) + ')';
  });

  this.link.attr('x1', function(d) { return d.source.x; })
    .attr('y1', function(d) { return d.source.y; })
    .attr('x2', function(d) { return d.target.x; })
    .attr('y2', function(d) { return d.target.y; });

};

CodeFlower.prototype.cleanup = function() {
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
};}

function getSize(d) {
  var bbox = this.getBBox(),
      cbbox = this.parentNode.getBBox(),
      scale = Math.min(cbbox.width/bbox.width, cbbox.height/bbox.height);
  d.scale = bbox.width;
};
