var page = require("webpage").create();

//var system = require('system');
//var url = "work/analysis/example_batch/Feature/demux_v/quality-plot.html"; //system.args[1];
//var img = "out.png";  //system.args[2];
//
var system = require('system'); //参数传递
var url = system.args[1]; //可以是在线网页地址, 也可以是本地html文件
var img = system.args[2]; //支持pdf/png等多种格式输出

//设置html视图窗口大小, 设置合适使得截图元素都在这个窗口显示
page.viewportSize = {
    width: 1600,
    height: 1200,
};

//设置截图区域
page.clipRect = {
    top: 20,
    left: 0,
    width: 800,
    height: 560
};

page.open(url, function(status) {
    if (status != "success" ){
        console.log("Unable to load the address!");
        phantom.exit();
    } else {
        //等待一会儿确保页面加载完成
        window.setTimeout(function(){
            //将页面保存为图片
            page.render(img);
            phantom.exit();
        }, 5000); //延迟时间可以根据实际情况调整
    }
});


