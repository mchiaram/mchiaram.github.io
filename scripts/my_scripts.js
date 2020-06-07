$(document).ready(function(){
  $("#menu-item").mouseover(function(){

    $("#menu-item").animate({font-size:30px},"fast");
   });
 $("#menu-item").mouseout(function(){
    $("#menu-item").animate({font-size:10px},"fast");
   });
})

$(document).ready(function(){
$('.multiple-items').slick({
 dots: true,
  infinite: true,
  speed: 700,
  slidesToShow: 1,
  centerMode: true,
  variableWidth: true,
  autoplay: true,
  autoplaySpeed: 1000,
	arrows: true
});
});


$(document).ready(function(){
$('.multiple-slides-intro').slick({
	dots: true,
  infinite: true,
  speed: 2000,
  slidesToShow: 1,
  centerMode: true,
  variableWidth: true,
  autoplay: true,
  autoplaySpeed: 2000,
	arrows: false
});
});

