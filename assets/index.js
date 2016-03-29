$(function(){
	$("address:contains('(a)')").each(function(){
		var $this = $(this), e = $this.text().replace(/ /g,'').replace('(a)','@');
		$this.empty().append($('<a>').attr('href','mailto:'+e).addClass('email').text(e));
	});
});