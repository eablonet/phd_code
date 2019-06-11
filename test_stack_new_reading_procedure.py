from packages.main import Stack as st

date = '20-11-2018'
serie = 5

s = st.Stack()
s.read_by_path(date, serie)
s.print_info()
s.current_image.show_image()
