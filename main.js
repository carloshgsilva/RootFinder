function fixedSize(x) {
    return x.toFixed(4);
}

const MAX_ITER = 16;
var PARAMS = {
    x0: 0,
    precision: 0,
    a: 0,
    b: 0,
    f: (x)=>x
}

function setParam(param, val){
    var f = undefined;
    if(param == 'f'){
        try{
            f = eval('(x)=>'+val);
            PARAMS.f = f;
            document.getElementById("field_f").classList.remove("failed");
        }catch(e){
            document.getElementById("field_f").classList.add("failed");
            return;
        }
    }else{
        PARAMS[param] = val;
    }

    saveParams();
    runAlgorithms();
}
function saveParams(){
    var p = {...PARAMS};
    p.f = p.f.toString();
    localStorage.params = JSON.stringify(p);
}
function loadParams(){
    if(localStorage.params == undefined){
        saveParams();
    }

    PARAMS = JSON.parse(localStorage.params);
    document.getElementById("field_f").value = PARAMS.f.substr(5); //removes '(x)=>'
    document.getElementById("field_x0").value = PARAMS.x0;
    document.getElementById("field_precision").value = PARAMS.precision;
    document.getElementById("field_a").value = PARAMS.a;
    document.getElementById("field_b").value = PARAMS.b;
    PARAMS.f = eval(PARAMS.f);
}

function isRoot(x, params){   
    var n = Math.floor(1.0/params.precision);
    var i = Math.trunc(Math.abs(params.f(x)*n));
    //console.log(`${Math.abs(i)} < ${1}`); //prec => 0.01 || n => 0.014111811976343258 || 
    return i < 1;
}

function bisectionMethod(params){
    var result = [];
    var a = params.a;
    var b = params.b;
    //Verificar se f(a) ou f(b) já é uma raíz aproximada
    if(isRoot(a, params)){
        return result;
    }
    if(isRoot(b, params)){
        return result;
    }
    if(params.f(a)*params.f(b) > 0){
        return result;
    }

    //Começar as iterações
    for(var i = 1; i < MAX_ITER; i++){
        //Calcular a bissecao x e o valor da funcao em x
        xi = (a + b) / 2;
        result.push({
            i: i,
            a: a,
            b: b,
            xi: xi,
            fa: params.f(a),
            fb: params.f(b),  
            fxi: params.f(xi),      
        });
        //console.log(`A: ${a} B: ${b} xi: ${xi} f(xi): ${params.f(xi)}`);
        //Verificar se x já é um zero: se sim abandonar iteracoes
        if (isRoot(xi, params)){
            break;
        }

        //Ajustar "a" ou "b" e o valor correspondente da função para que o inter-
        //Intervalo [a,b] diminua mas continue contendo a raiz
        if (params.f(params.a)*params.f(xi) < 0){
            b = xi;
        }
        else{
            a = xi;
        }       
    }
    return result;
}

const EPSILON = 0.00001
function dx(x, params){
    var y0 = params.f(x-+EPSILON);
    var y1 = params.f(x+EPSILON);
    if(y1 == y0)return 1.0;
    return (y1-y0)/(2.0*EPSILON);
}

function newtonRaphson(params){
    var result = [];
    var last_k = params.x0;

    for(var i = 0; i < MAX_ITER; i++){
        var fx_val = params.f(last_k, params.precision);
        var dx_val = dx(last_k, params);
        var k = last_k - fx_val/dx_val;
        if(isNaN(k))break;
        result.push({
            i: i+1,
            last_k: last_k,
            fx_val: fx_val,
            dx_val: dx_val,
            k: k,
            f_k: params.f(k),         
        });
        if(Math.abs(params.f(k)) < params.precision){
            break;
        }
        last_k = k;
    }
    return result;
}

function renderResults(newtonRaphson, bisec){
    var result = document.getElementById("result_table");
    var result2 = document.getElementById("result_table2");
    var solution = document.getElementById("result_solution");
    var solution2 = document.getElementById("result_solution2");

    result.innerHTML = `    
        <tr class="table_header">
            <th>$$ i $$</th>
            <th>$$ x_i $$</th>
            <th>$$ f(x_i) $$</th>
            <th>$$ f'(x_i) $$</th>
            <th>$$ x_{i + 1} $$</th>
            <th>$$ f(x_{i + 1}) $$</th>
        </tr>
        ${newtonRaphson.map(r => `
            <tr>
                <td>${r.i}</td>
                <td>${fixedSize(r.last_k)}</td>
                <td>${fixedSize(r.fx_val)}</td>
                <td>${fixedSize(r.dx_val)}</td>
                <td>${fixedSize(r.k)}</td>
                <td>${fixedSize(r.f_k)}</td>
            </tr>            
        `).join("")}
    `;

    result2.innerHTML = `    
        <tr class="table_header">
            <th>$$ i $$</th>
            <th>$$ a $$</th>
            <th>$$ b $$</th>
            <th>$$ x_i $$</th>
            <th>$$ f(a) $$</th>
            <th>$$ f(b) $$</th>
            <th>$$ f(x_i) $$</th>
        </tr>
        ${bisec.map(r => `
        <tr>
            <td>${r.i}</td>
            <td>${fixedSize(r.a)}</td>
            <td>${fixedSize(r.b)}</td>
            <td>${fixedSize(r.xi)}</td>
            <td>${fixedSize(r.fa)}</td>
            <td>${fixedSize(r.fb)}</td>
            <td>${fixedSize(r.fxi)}</td>
        </tr>            
        `).join("")}
    `;

    document.getElementById("result_solution").innerHTML = `$$${newtonRaphson.map(r => `
         X_{${r.i}} = ${fixedSize(r.last_k)} - \\frac{${fixedSize(r.fx_val)}}{${fixedSize(r.dx_val)}} \\\\ \\ \\\\
    `) .join("")} \\\\ x = ${newtonRaphson[newtonRaphson.length - 1].k}$$`;

    if(typeof renderMathInElement == "function"){
        renderMathInElement(solution);
        renderMathInElement(solution2);
        renderMathInElement(result);
        renderMathInElement(result2);
    }

    //Find bounding box
    var bbox = [999999,999999,-999999,-999999];
    for(var r of newtonRaphson){
        bbox[0] = Math.min(bbox[0], r.last_k);
        bbox[1] = Math.min(bbox[1], r.fx_val);
        bbox[2] = Math.max(bbox[2], r.last_k);
        bbox[3] = Math.max(bbox[3], r.fx_val);
    }
    var bbox_eps = Math.max(Math.abs(bbox[0]-bbox[2]), Math.abs(bbox[1]-bbox[3]))*0.1;

    //init board
    const board = JXG.JSXGraph.initBoard('graph', { 
        boundingbox: [bbox[0]-bbox_eps, bbox[3]+bbox_eps, bbox[2]+bbox_eps, bbox[1]-bbox_eps], axis:true, 
        keepAspectRatio: false,
    });

    //Create points
    var points = [];
    for(var r of newtonRaphson){
        points.push(board.create('point', [''+r.last_k,''+r.fx_val], {style:6, name:''+r.i}));
    }
    for(var i=0;i<points.length-1;i++){
        var a = points[i];
        var b = points[i+1];
        board.create('arrow', [a, b], {color: "#000a", dash: 1});
    }

    // Macro function plotter
    function addCurve(board, func, atts) {
        var f = board.create('functiongraph', [func], atts);
        return f;
    }
    // Simplified plotting of function
    function plot(func, atts) {
        if (atts==null) { return addCurve(board, func, {strokewidth:2});
        } else { return addCurve(board, func, atts); }    
    }

    function func(x) {
        return PARAMS.f(x);
     }
     c = plot(func);
     // Derivative:
     var g = JXG.Math.Numerics.D(PARAMS.f);
     plot(g,{strokecolor:'black'});
     
}

function runAlgorithms(){
    var newtonRaphsonResult = newtonRaphson(PARAMS);
    var bisecResult = bisectionMethod(PARAMS);
    renderResults(newtonRaphsonResult, bisecResult);
}




function expectNR(params, expected){
    var result = newtonRaphson(params);
    var k = result[result.length-1].k;
    if(Math.abs(k-expected) <= params.precision){
        console.info(`Success! x0: ${params.x0} ${params.f} -> ${k} == ${expected}`);
    }else{
        console.error(`Failed! x0: ${params.x0} ${params.f} -> ${k} != ${expected}`);
        console.warn(result);
    }
}
function expectBisection(params, expected){
    var result = bisectionMethod(params);
    var xi = result [result.length-1].xi;
    if(Math.abs(xi-expected) <= params.precision){
        console.info(`Success! a: ${params.a} b: ${params.b} f: ${params.f}`, "\n", ` result: ${xi} == exptected: ${expected}`);
    }else{
        console.error(`Failed! a: ${params.a} b: ${params.b} f: ${params.f}`, "\n" , ` result: ${xi} != exptected: ${expected}`);
        console.warn(result);
    }
}
function runTests(){
    //f(x) = x^3 + 2x - 5
    expectBisection({a: 1, b: 2, f: (x)=>x*x*x + 2*x - 5, precision: 0.01}, 1.3281);
    expectNR({x0: 1, f: (x)=>x*x*x + 2*x - 5, precision: 0.001}, 1.3283);

    //f(x) = x^2 - 4 + e^x
    expectBisection({a: -2, b: -1, f: (x)=>x*x - 4 + Math.pow(Math.E, x), precision: 0.01}, -1.9648);
    expectBisection({a: 1, b: 2, f: (x)=>x*x - 4 + Math.pow(Math.E, x), precision: 0.01}, 1.0586);
    expectNR({x0: -2, f: (x)=>x*x - 4 + Math.pow(Math.E, x), precision: 0.001}, -1.9646);
    expectNR({x0: 1, f: (x)=>x*x - 4 + Math.pow(Math.E, x), precision: 0.001}, 1.0580);
    
    //f(x) = sqrt(x+9) - (x^2)/2
    expectBisection({a: -3, b: -2, f: (x)=>Math.sqrt(x+9) - (x*x)/2, precision: 0.01}, -2.2812);
    expectBisection({a: 2, b: 3, f: (x)=>Math.sqrt(x+9) - (x*x)/2, precision: 0.01}, 2.6094);
    expectNR({x0: -3, f: (x)=>Math.sqrt(x+9) - (x*x)/2, precision: 0.001}, -2.2772);
    expectNR({x0: 2, f: (x)=>Math.sqrt(x+9) - (x*x)/2, precision: 0.001}, 2.6105);

    //f(x) = 5cos(x) - x^2
    expectBisection({a: -2, b: -1, f: (x)=>5*Math.cos(x) - x*x, precision: 0.01}, -1.2520);
    expectBisection({a: 1, b: 2, f: (x)=>5*Math.cos(x) - x*x, precision: 0.01}, 1.2520);
    expectNR({x0: -2, f: (x)=>5*Math.cos(x) - x*x, precision: 0.001}, -1.2519);
    expectNR({x0: 1, f: (x)=>5*Math.cos(x) - x*x, precision: 0.001}, 1.2521);

    //f(x) = 4sen(x) - e^x
    expectBisection({a: 1, b: 2, f: (x)=>4*Math.sin(x) - Math.pow(Math.E, x), precision: 0.01}, 1.3672);
    expectNR({x0: 1, f: (x)=>4*Math.sin(x) - Math.pow(Math.E, x), precision: 0.001}, 1.3650);
}

//runTests();

document.body.onload=()=>{
    loadParams();
    runAlgorithms();
}
