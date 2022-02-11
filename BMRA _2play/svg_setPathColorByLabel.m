% xDoc = svg_setPathColorByLabel(xDoc, label, color)
%   Set color of inkscape svg path by its inkscape label
%   xDoc ... DOM node of the svg file
%   label ... string with label name
%   color ... 1x3 vector, rbg color

function [xDoc success] = svg_setPathColorByLabel(xDoc, label, color)

allPathitems = xDoc.getElementsByTagName('path');

success = false;

for i= 0:allPathitems.getLength-1
    thisPathitem = allPathitems.item(i);
    thisAttr = thisPathitem.getAttributeNode("inkscape:label");
    if ~isempty(thisAttr)
        thisLabel = string(thisAttr.getValue);
        if thisLabel==join(["#" label], '')
            thisStyle = thisPathitem.getAttributeNode("style");
            thisString = string(thisStyle.getValue);
            %also: "stroke-opacity:1"
            before = extractBefore(thisString,"stroke:");
            after = extractAfter(thisString,"stroke:");
            after = extractAfter(after,";");
            myStrokecolor = rgb2hex(color);
            newString = join(["stroke:" myStrokecolor ";" before after],'');
            thisStyle.setValue(newString)
            %sprintf('%d yes',i)
            success = true;
        else
            %sprintf('%d no',i)
        end
    end
end

if ~success
    fprintf('  %s not found in svg\n',label)
end