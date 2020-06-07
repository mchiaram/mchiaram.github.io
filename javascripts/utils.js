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




