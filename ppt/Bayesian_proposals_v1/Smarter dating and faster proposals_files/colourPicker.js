


function getColourPicker(element, id, currentCol = null, resolve = function(col) { }) {

	var html = getColourPickerHTML(id, element.offset().left, element.offset().top, currentCol, resolve);
	$("body").append(html);

	



}



function pickColour(id, resolve, colEle) {
	$("#colourpicker_" + id).remove();
	var col = $(colEle).css('background-color');
	resolve(col, id);
}




function getColourPickerHTML(id, x, y, currentCol = "black", resolve = function(col, id) { }) {


	var colourPickHtml = `
	<div id="colourpicker_` + id + `" style="top:` + y + `; left:` + x + `" class="colourpickerrow">
		<td colspan=3>
			<ul class="flex-container thicklines" style="font-size:100%;">`;			
			
	
	var colMatch = false;
	for (var colNum = 0; colNum < DEFAULT_COLOURS.length; colNum++){
		var col = DEFAULT_COLOURS[colNum].toUpperCase();
		//console.log("col", col, currentCol);
		
		var addClass = "";
		if (col == currentCol){
			addClass = "selected";
			colMatch = true;
		}

		colourPickHtml += `
				<li class="flex-item">
					<span class="colourbox ` + addClass + `" title="` + col + `" onclick="pickColour('` + id + `',` + resolve + `, this);" style="background-color:` + col + `"></span>
				</li>`;
		
	}		
	
	// Custom hex codes
	/*
	colourPickHtml += `
				<li class="flex-item" style="margin-left:10px" title="Enter hex code">
					Hex code: <input id="geneColHexCode` + id + `" class="numberinput" onchange="pickColour('` + id + `',` + resolve + `, this);" style="width:6em" value="#">
				</li>`;
	*/
	colourPickHtml += `</ul></span>`;

	return colourPickHtml;
	

}









