var cityDropDown = document.getElementById('facilityCity');
var tankTypeDropDown = document.getElementById('tankType');
var tankConditionDropDown = document.getElementById('tankCondition');
var tankColorDropDown = document.getElementById('tankColor');
var tankShadeDropDown = document.getElementById('tankShade');
var tankLiquidDropDown = document.getElementById('tankLiquid');

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