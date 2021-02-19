roses2020 unit01 readme

environment:
- Use the default 'roses' environment.

notebooks:
- The "follow along" notebook is completed with the lecture video.
- The "full" notebook is a solution and should execute without exceptions.
- the "lab" notebook contains prompts for a student to complete.

data:
- 4 miniseed files are included.

faq/notes:
- Refer to the official Obspy tutorial for more info: docs.obspy.org/tutorial/
- Consider using Python's f-strings when building a string from variables (like nslc codes).
- Obspy's 'get_stations' and 'attach_response' can be used if 'get_waveforms(...,attach_response=True) does not work as expected.
- Obspy's 'read_inventory' function is used to attach response from a local file.

