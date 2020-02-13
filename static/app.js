var results="";
var cityDropDown = document.getElementById('facilityCity');
var tankTypeDropDown = document.getElementById('tankType');
var tankConditionDropDown = document.getElementById('tankCondition');
var tankColorDropDown = document.getElementById('tankColor');
var tankShadeDropDown = document.getElementById('tankShade');
var contents = document.getElementsByName('contents');
var FacilityName = document.getElementById('FacilityName');
var FacilityAddress = document.getElementById('FacilityAddress');
var facilityCity = document.getElementById('facilityCity');
var TankName = document.getElementById('TankName');
var TankLongitude = document.getElementById('TankLongitude');
var TankLatitude = document.getElementById('TankLatitude');
var tankLiquidDropDown = document.getElementById('tankLiquid');
var TankHeight = document.getElementById('TankHeight');
var TankDiameter = document.getElementById('TankDiameter');
var AvgLiqHeight = document.getElementById('AvgLiqHeight');
var AnnThroughput = document.getElementById('AnnThroughput');
var AvgLowerLiquid = document.getElementById('AvgLowerLiquid');
var AvgHigherLiquid = document.getElementById('AvgHigherLiquid');
var tankType = document.getElementById('tankType');
var tankLiquid = document.getElementById('tankLiquid');
var tankCondition = document.getElementById('tankCondition');
var tankColor = document.getElementById('tankColor');
var tankShade = document.getElementById('tankShade');
var content = [];
var quantity = [];
for (var i = 0; i < input_chem_array.length; i++) {
    var unitOption = input_chem_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;

    for (var j = 0; j < contents.length; j++) {

        contents[j].appendChild(element);
}
}


for (var i = 0; i < input_met_array.length; i++) {
    var unitOption = input_met_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    cityDropDown.appendChild(element);
}

for (var i = 0; i < input_condition_color_array.length; i++) {
    var unitOption = input_condition_color_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    tankColorDropDown.appendChild(element);
}

for (var i = 0; i < input_condition_shade_array.length; i++) {
    var unitOption = input_condition_shade_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    tankShadeDropDown.appendChild(element);
}

for (var i = 0; i < input_condition_condition_array.length; i++) {
    var unitOption = input_condition_condition_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    tankConditionDropDown.appendChild(element);
}

for (var i = 0; i < input_tank_type_array.length; i++) {
    var unitOption = input_tank_type_array[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    tankTypeDropDown.appendChild(element);
}

for (var i = 0; i < input_contents_type.length; i++) {
    var unitOption = input_contents_type[i];
    var element = document.createElement("option");
    element.textContent = unitOption;
    element.value = unitOption;
    tankLiquidDropDown.appendChild(element);
}


$(document).ready(function () {

  


    $("#add_row").on("click", function () {
        // Dynamic Rows Code

        // Get max row id and set new id
        var newid = 0;
        $.each($("#tab_logic tr"), function () {
            if (parseInt($(this).data("id")) > newid) {
                newid = parseInt($(this).data("id"));
            }
        });
        newid++;

        var tr = $("<tr></tr>", {
            id: "addr" + newid,
            "data-id": newid
        });

        // loop through each td and create new elements with name of newid
        $.each($("#tab_logic tbody tr:nth(0) td"), function () {
            var td;
            var cur_td = $(this);

            var children = cur_td.children();

            // add new td and element if it has a nane
            if ($(this).data("name") !== undefined) {
                td = $("<td></td>", {
                    "data-name": $(cur_td).data("name")
                });

                var c = $(cur_td).find($(children[0]).prop('tagName')).clone().val("");
                c.attr("name", $(cur_td).data("name") + newid);
                c.appendTo($(td));
                td.appendTo($(tr));
            } else {
                td = $("<td></td>", {
                    'text': $('#tab_logic tr').length
                }).appendTo($(tr));
            }
        });

       
        // add the new row
        $(tr).appendTo($('#tab_logic'));

        $(tr).find("td button.row-remove").on("click", function () {
            $(this).closest("tr").remove();
        });
    });




    // Sortable Code
    var fixHelperModified = function (e, tr) {
        var $originals = tr.children();
        var $helper = tr.clone();

        $helper.children().each(function (index) {
            $(this).width($originals.eq(index).width())
        });

        return $helper;
    };


    $("#add_row").trigger("click");
    

   
   $('#myModal').on('shown.bs.modal', function () {
   
    $('#myInput').trigger('focus')
  });
   function getdata()
   {

        getquntityandcontents();
        var POST_JSON = {};

        POST_JSON["input_city"] = facilityCity.value;
        POST_JSON["input_tank"] = [parseFloat(TankHeight.value), parseFloat(AvgLiqHeight.value), parseFloat(TankDiameter.value)];
        POST_JSON["input_contents"] = [parseFloat(AnnThroughput.value), tankLiquid.value, parseFloat(AvgLowerLiquid.value), parseFloat(AvgHigherLiquid.value)];
        POST_JSON["input_chem"] = content;
        POST_JSON["input_qty"] = quantity;
        POST_JSON["input_default"] = [0.0625, 1];
        POST_JSON["input_condition"] = [tankColor.value, tankShade.value, tankCondition.value];
        POST_JSON["input_tank_type"] = tankType.value;
        POST_JSON["input_tank_name"] = TankName.value;
        POST_JSON["input_facility_name"] = FacilityName.value;
        POST_JSON["input_facility_address"] = FacilityAddress.value;
        POST_JSON["input_tank_geo"] = [parseFloat(TankLongitude.value), parseFloat(TankLatitude.value)];
        // Test 
       // POST_JSON =   { "input_city": "Los Angeles AP, California", "input_tank": [12, 9, 8], "input_contents": [8450, "other stocks", 11.5, 4.5], "input_chem": ["Cyclohexane", "Benzene", "Toluene"], "input_qty": [101, 2812, 258], "input_default": [0.0625, 1], "input_condition": ["White", "None", "Average"], "input_tank_type": "Vertical", "input_tank_name": "Tank 100019", "input_facility_name": "Tank 100019", "input_facility_address": "2101 Pacific Coast Hwy, Wilmington, CA 90744", "input_tank_geo": [33.793548, -118.23332] }
        
        $.ajax
            ({
                type: "POST",
                //the url where you want 
                url: '/api/vfrtk',
                dataType: 'json',
                async: false,
                contentType: 'application/json;charset=UTF-8',
                data:JSON.stringify(POST_JSON, null, '\t'),
                success: function (data) {

                    results =data;

                    CreateDetailsTable();
                    CreateSummaryTable();
                }
            })
         
         
        console.log(JSON.stringify(POST_JSON));
        quantity = [];
        content = [];
      
  //  });
        }
      
    
    function getquntityandcontents() {

        $("#tab_logic").find('tbody > tr').each(function (i) {

            if (i == 0) return;
            var $fieldset = $(this);
            qty = $('input:text:eq(0)', $fieldset).val();
            cont = $('select option:selected:eq(0)', $fieldset).val();

            quantity.push(parseFloat(qty));
            content.push(cont);
        });

        }
    

        function CreateDetailsTable(){
            var data =  results.detail; 
            var table = document.createElement("table");
            table.setAttribute('class', 'table table-bordered');
            var colorder=
                {
                 "antoine_coef_a": 8,
                  "antoine_coef_b": 9,
                  "antoine_coef_c": 10,
                  "cas_no": 5,
                  "comp_amt": 7,
                  "comp_mole": 14,
                  "comp_mole_xi":15,
                  "comp_partial": 16,
                  "comp_partial_tln":19,
                  "comp_partial_tlx": 20,
                  "comp_vap_mole_frac": 17,
                  "comp_vapor_mw_xi": 18,
                  "comp_vp":11,
                  "comp_vp_tln": 13,
                  "comp_vp_tlx": 12,
                  "component": 4,
                  "facility_name": 1,
                  "mw": 6,
                  "stand_loss_xi": 23,
                  "tank_name": 2,
                  "tank_type": 3,
                  "total_loss_xi": 25,
                  "vap_mole_xi": 21,
                  "vap_wt_xi":22,
                  "work_loss_xi":24
      
            }
          
           results.order = colorder;
           var columns= results.order;
            var col = [];
           for(var k in columns)
           {

            col[columns[k]-1] = k;
          

           }
                 for (var i =0; i< col.length ; i++)
                 {
                    tr = table.insertRow(-1); 
                    var tabCell = tr.insertCell(-1);
                    tabCell.innerHTML = col[i];

                    for (var prop in data) {                  
                   
                    var tabCell2 = tr.insertCell(-1);
                    tabCell2.innerHTML = data[prop][col[i]];


                    }


                 }             

      var divContainer = document.getElementById("emissiondetail");
        divContainer.innerHTML = "";
        divContainer.appendChild(table);

        }
        function CreateSummaryTable()
        {

            var data =  results.summary; 
            var table = document.createElement("table");
            table.setAttribute('class', 'table table-bordered');
            var col = [];
            var transposecol=[];
            for (var key in data) {
                    if (col.indexOf(key) === -1) {
                        col.push(key);
                    }

                    for(var k in data[key])
                    {

                        if (transposecol.indexOf(k) === -1) {
                            transposecol.push(k);
                        }
                    }
                
            }


            tr = table.insertRow(-1); 
            var tabCell = tr.insertCell(-1);
            tabCell.innerHTML = "";
            for (var i =0; i< col.length ; i++)
            {
               
               var tabCell = tr.insertCell(-1);
               tabCell.innerHTML = col[i];
            }
        
        
        
        for (var j =0; j< transposecol.length ; j++)
        {   
            tr = table.insertRow(-1); 
            var tabCell = tr.insertCell(-1);
            tabCell.innerHTML =transposecol[j];  
        for (var i =0; i< col.length ; i++)
        {         
               
                var tabCell = tr.insertCell(-1);
                tabCell.innerHTML = data[col[i]][transposecol[j]];    

                
                
        }
        }


        var divContainer = document.getElementById("summary");
        divContainer.innerHTML = "";
        divContainer.appendChild(table);
    }
        
    

 var modal = document.getElementById("myModal");

 // Get the button that opens the modal
 var btn = document.getElementById("submit");

 // Get the <span> element that closes the modal
 var span = document.getElementsByClassName("close")[0];

 // When the user clicks on the button, open the modal
 btn.onclick = function() {
  getdata();
  modal.style.display = "block";
 }

 // When the user clicks on <span> (x), close the modal
 span.onclick = function() {
  modal.style.display = "none";
 }

// When the user clicks anywhere outside of the modal, close it
 window.onclick = function(event) {
  if (event.target == modal) {
    modal.style.display = "none";
  }
}

});



