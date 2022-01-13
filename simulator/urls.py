from django.conf.urls import url
from simulator.views import home, simulator, default

urlpatterns = [
	url(r'^$',home),
	url(r'index',home),
	url(r'simulator', simulator),
	url(r'default', default)
]