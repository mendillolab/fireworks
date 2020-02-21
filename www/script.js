var nodes = null;
var nodes_import = null;
var nodes_array = [];

var edges = null;
var edges_import = null;
var edges_array = [];

var network = null;
var options = null;

// get node data from R
Shiny.addCustomMessageHandler('svg_handler_nodes', function(message){
  nodes_import = Object.entries(message);
  nodes_array = [];
  nodes_import.forEach(([key, value]) => {
    nodes_array.push(value)
  })
  console.log(nodes_array)
});

// get edge data from R
Shiny.addCustomMessageHandler('svg_handler_edges', function(message){
  edges_import = Object.entries(message);
  edges_array = [];
  edges_import.forEach(([key, value]) => {
    edges_array.push(value)
  })
  console.log(edges_array)
});

// update & re-draw network
Shiny.addCustomMessageHandler('svg_handler_update', function(message){
  draw();
  exportSvg()
});

function draw() {

  nodes = new vis.DataSet(nodes_array);
  edges = new vis.DataSet(edges_array);
  // // create an array with nodes
  // nodes = new vis.DataSet([
  //   { id: 1, label: "C16orf72", x: 15, y: 15},
  //   { id: 2, label: "HUWE1",    x: 10,y: 25}
  // ]);
  //
  // // create an array with edges
  // edges = new vis.DataSet([
  //   { from: 1, to: 2 },
  // ]);

  // create a network
  var container = document.getElementById("network_svg_canvas");
  var data = {
    nodes: nodes,
    edges: edges
  };
  options = {
    height: '600px',
    width: '800px',
    physics: false,
    interaction: {
      dragNodes: false,
      zoomView: false,
      dragView: false
    },
    nodes: {
      shape: 'dot',
      font: {
        vadjust:-50,
        size: 30,
        strokeWidth:0,
        strokeColor:"#000000"
      }
    },
    groups: {
      source: {
        size: 45,
        font: {
          size:45,
          vadjust: -85
        }
      }
    },
    edges: {
      smooth: {
        enabled: true,
        type: "dynamic",
        roundness: 1.0
      },
      length:15
    },
  };
  network = new vis.Network(container, data, options);


  console.log(network)
}

// canvas2svg conversion code
C2S.prototype.circle = CanvasRenderingContext2D.prototype.circle;
C2S.prototype.square = CanvasRenderingContext2D.prototype.square;
C2S.prototype.triangle = CanvasRenderingContext2D.prototype.triangle;
C2S.prototype.triangleDown = CanvasRenderingContext2D.prototype.triangleDown;
C2S.prototype.star = CanvasRenderingContext2D.prototype.star;
C2S.prototype.diamond = CanvasRenderingContext2D.prototype.diamond;
C2S.prototype.roundRect = CanvasRenderingContext2D.prototype.roundRect;
C2S.prototype.ellipse_vis = CanvasRenderingContext2D.prototype.ellipse_vis;
C2S.prototype.database = CanvasRenderingContext2D.prototype.database;
C2S.prototype.arrowEndpoint = CanvasRenderingContext2D.prototype.arrowEndpoint;
C2S.prototype.circleEndpoint = CanvasRenderingContext2D.prototype.circleEndpoint;
C2S.prototype.dashedLine = CanvasRenderingContext2D.prototype.dashedLine;

function exportSvg()
{
    var networkContainer = network.body.container;
    var ctx = new C2S({width: networkContainer.clientWidth, height: networkContainer.clientWidth, embedImages: true});

    var canvasProto = network.canvas.__proto__;
    var currentGetContext = canvasProto.getContext;
    canvasProto.getContext = function()
    {
        return ctx;
    }
    var svgOptions = {
        nodes: {
            shapeProperties: {
                interpolation: false //so images are not scaled svg will get full image
            },
            scaling: { label: { drawThreshold : 0} },
            font:{color:'#000000'}
        },
        edges: {
            scaling: { label: { drawThreshold : 0} }
        }
    };
    network.setOptions(svgOptions);
    network.redraw();
    network.setOptions(options);
    canvasProto.getContext = currentGetContext;
    ctx.waitForComplete(function()
        {
            var svg = ctx.getSerializedSvg();
            showSvg(svg);

        });
}

function showSvg(svg)
{


    var svgBlob = new Blob([svg], {type: 'image/svg+xml'});
    openBlob(svgBlob, "network.svg");
    //svgElementToPdf(svg, pdf, {
    //  scale: 72/96, // this is the ratio of px to pt units
    //  removeInvalid: true // this removes elements that could not be translated to pdf from the source svg
    //});
    // var pdfBlob = new Blob([pdf], {type: 'image/pdf'});
    // openBlob(pdfBlob, "network.pdf");
}

function openBlob(blob, fileName)
{
if(window.navigator && window.navigator.msSaveOrOpenBlob)
    {

        //blobToDataURL(blob, function(dataurl){window.open(dataurl);});
        window.navigator.msSaveOrOpenBlob(blob,fileName);
    }
    else
    {
 var a = document.getElementById("blobLink");
 if(!a)
 {
   a = document.createElement("a");
   document.body.appendChild(a);
   a.setAttribute("id", "blobLink");
   a.style = "display: none";
 }
 var data = window.URL.createObjectURL(blob);
 a.href = data;
 a.download = fileName;
 a.click();
 setTimeout(function()
   {
   // For Firefox it is necessary to delay revoking the ObjectURL
   window.URL.revokeObjectURL(data);
   }
   , 100);
    }
}
