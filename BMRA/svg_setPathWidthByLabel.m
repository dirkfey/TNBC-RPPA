% xDoc = svg_setPathWidthByLabel(xDoc, label, width)
%   Set color of inkscape svg path by its inkscape label
%   xDoc ... DOM node of the svg file
%   label ... string with label name
%   width ... double (in px)

function [xDoc success] = svg_setPathWidthByLabel(xDoc, label, width)

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
            before = extractBefore(thisString,"stroke-width:");
            after = extractAfter(thisString,"stroke-width:");
            after = extractAfter(after,";");
            newString = join(["stroke-width:" width "px;" before after],'');
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