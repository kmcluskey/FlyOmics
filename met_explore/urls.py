from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='met_explore_index'),
    path('metabolite_list', views.MetaboliteListView.as_view(), name='metabolite_list'),
]

