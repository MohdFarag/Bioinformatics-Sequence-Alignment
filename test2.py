import flet
from flet import *
from urllib.parse import urlparse


def main(page: Page):
	page.title = "routing app"
	youparams = "watermelon"

	def route_change(route):
		# CLEAR ALL PAGE
		page.views.clear()
		page.views.append(
			View(
				"mm/jjlllllllj",
				[
				AppBar(title=Text("Home Page", size=30,
					color="white"
					),
					bgcolor="blue",
					),
					
                    # PAGE ROUTE IS PATH YOU URL HERE
					Text(page.route),
					ElevatedButton(
						"Go to Second Page",
						on_click=lambda _: page.go(f"/secondpage/{youparams}")
						)
				])
			)

		# GET param from home page
		param = page.route
		# THIS IS GET VALUE AFTER /secondpage/THIS RES HERE
		res = urlparse(param).path.split("/")[-1]
		print(f"test res is : {res}")

		if page.route == f"/secondpage/{res}":
			page.views.append(
				View(
				# IF URL ACCESS HERE THEN PUSH TO VIEW HERE
				f"/secondpage/{res}",
				[
					# LIKE BEFORE
					# RENDER YOU PAGE HERE
					# PAGE ROUTE IS PATH YOU URL HERE
					Text(page.route),
					Text(f"you params is {res}"),
					ElevatedButton(
						"BACK TO HOME PAGE",
						on_click=lambda _: page.go("/")
						)
                        ]))
	page.update()

	def view_pop(view):
		page.views.pop()
		myview = page.views[-1]
		page.go(myview.route)

	page.on_route_change = route_change
	page.on_view_pop = view_pop
	page.go(page.route)
	p = TemplateRoute(page.route)
	if p.match("/second/:id"):
		print("you here ", p.id)
	else:
		print("whatever")

flet.app(target=main)
