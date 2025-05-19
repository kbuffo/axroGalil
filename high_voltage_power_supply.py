"""
Information and specs of the high voltage power supply.
"""

# Info generated when running Keysight Connection Expert of Keysight IO Library Suite
psu_model = 'N5771A'
psu_manufacturer = 'Agilent Technologies'
psu_serial_number = 'US15H8632R'
psu_firmware_version = 'B.00.08.309'
psu_web_link = 'https://www.keysight.com/de/de/product/N5771A/power-supply-300v-5a-1500w.html'
psu_visa_address = 'USB0::0x0957::0xA807::US15H8632R::0::INSTR'
psu_alias = 'USBInstrument1'
psu_sicl_address = 'usb0[2391::43015::US15H8632R::0]'

default_psu_min_voltage = 0. # V, the default minimum output voltage allowed
default_psu_max_voltage = 314.2 # V, the default maximum output voltage allowed
psu_min_voltage = 0. # V, the software constrained minimum output voltage allowed
psu_max_voltage = 310. # V, the software constrained maximum output voltage allowed

default_psu_min_current = 0. # A, the default minimum output current allowed
default_psu_max_current = 5. # A, the default maximum output current allowed
psu_min_current = 0. # A, the software constrained minimum output current allowed
psu_max_current = 0.04 # A, the software constrained maximum output current allowed