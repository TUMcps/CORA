const dur = 10;
const svgAnms = document.getElementsByClassName("svg-animate");
for (var i=0; i<svgAnms.length; i++) {
    var children = svgAnms[i].children;
    const delay = dur/children.length/5;
    for (let j = 0; j < children.length; j++) {
        children[j].innerHTML = "<animate attributeName=\"opacity\" values=\"0;0.5;1;1;1;1;0;0;0;0\" dur=\"" + dur + "s\"  begin=\"" + delay*j + "s;op.end+" + delay*j + "s\" repeatCount=\"indefinite\" />";
    }
}